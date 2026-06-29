#!/usr/bin/env python3
r"""Regularized finite adjoint-Green kernel model.

The raw pseudoinverse coefficients in ``adjoint_green_boundary_diagnostic.py``
prove the finite active range identity

    E_active = C_N R_N,

but they are only grid coefficients.  Differentiating them directly is not a
valid approximation to the adjoint Green boundary-value problem.

This script constructs smoother finite kernels instead.  For every active
source row it finds a polynomial density K(a) on I=[a_-,a_+] which satisfies

    int_I K(a) Lambda_a(v_j) da = E_active(v_j)

on the two-dimensional active source space, while minimizing a fixed Sobolev-
type polynomial norm.  It then tests the resulting Green identity

    int_I K Lambda(f)
      = int_I f P^*K + [B_P[K,f]]_{a_-}^{a_+}

on the same active space.  This is still finite evidence, not the continuum
proof, but it is the correct regularized object for a compactness/convergence
argument.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_boundary_diagnostic import (  # noqa: E402
    active_indices,
    boundary_poly_row,
    density_derivatives,
    integral_poly_row,
    row_norm,
    trace_args,
    weights_for_centers,
)
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_active_range_inclusion import inverse_spd  # noqa: E402
from global_trace_observability_gap import (  # noqa: E402
    child_args,
    f,
    fmt,
    frob_norm,
    source_nodes,
    stack_source_rows,
)
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from lambda_differential_closure import coefficient_derivatives  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    poly_derivative_value,
    poly_value,
    shifted_legendre_polys,
    trace_matrix,
)
from source_concomitant_membership import sample_trace_field  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def trace_matrix_for_samples(args, polys, samples: int):
    targs = SimpleNamespace(**vars(args))
    targs.global_constraints = samples
    targs.constraints = samples
    targs.constraint_min = args.constraint_min
    targs.constraint_max = args.constraint_max
    return trace_matrix(polys, make_qargs(targs))


def field_derivatives(args, centers, max_q: int):
    field_args = SimpleNamespace(**vars(args))
    field_args.s_min = args.constraint_min
    field_args.s_max = args.constraint_max
    field_args.samples = len(centers)
    centers_field, vecs, min_gap, min_cos = sample_trace_field(field_args)
    max_center_error = max(abs(centers[i] - centers_field[i]) for i in range(len(centers)))
    if args.stencil_width <= max_q:
        raise SystemExit("--stencil-width must exceed max derivative order")
    return (
        [
            coefficient_derivatives(centers, vecs, idx, max_q, args.stencil_width)
            for idx in range(len(centers))
        ],
        min_gap,
        min_cos,
        max_center_error,
    )


def polynomial_design(kpolys, centers, weights, trace_active, lo):
    active_dim = trace_active.cols
    design = mp.matrix(active_dim, len(kpolys))
    for col, poly in enumerate(kpolys):
        for mode in range(active_dim):
            design[mode, col] = mp.fsum(
                weights[i] * poly_value(poly, centers[i] - lo) * trace_active[i, mode]
                for i in range(len(centers))
            )
    return design


def regularized_coefficients(design, target, degree: int, smooth_order: int, smooth_weight):
    """Minimum H-norm exact solution of design*c=target."""
    h_inv = [
        1
        / (
            mp.mpf("1")
            + smooth_weight
            * (mp.mpf(j) / max(1, degree - 1)) ** (2 * smooth_order)
        )
        for j in range(degree)
    ]
    active_dim = design.rows
    reduced = mp.matrix(active_dim)
    for i in range(active_dim):
        for j in range(active_dim):
            reduced[i, j] = mp.fsum(
                design[i, m] * h_inv[m] * design[j, m] for m in range(degree)
            )
    rhs = mp.matrix(active_dim, 1)
    for i in range(active_dim):
        rhs[i] = target[i]
    alpha = mp.lu_solve(reduced, rhs)
    coeffs = [
        h_inv[m] * mp.fsum(design[i, m] * alpha[i] for i in range(active_dim))
        for m in range(degree)
    ]
    return coeffs


def kernel_derivatives(kpolys, coeffs, point, lo, max_q):
    y = point - lo
    return [
        mp.fsum(coeffs[m] * poly_derivative_value(kpolys[m], y, q) for m in range(len(kpolys)))
        for q in range(max_q + 1)
    ]


def kernel_values(kpolys, coeffs, centers, lo):
    return [
        mp.fsum(coeffs[m] * poly_value(kpolys[m], center - lo) for m in range(len(kpolys)))
        for center in centers
    ]


def list_norm(values):
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=11)
    parser.add_argument("--active-tol", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--trace-samples", type=int, default=41)
    parser.add_argument("--kernel-degree", type=int, default=12)
    parser.add_argument("--smooth-order", type=int, default=4)
    parser.add_argument("--smooth-weight", default="1e-6")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--stencil-width", type=int, default=13)
    parser.add_argument("--json-out", default="adjoint_green_bvp_regularized.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    active_tol = mp.mpf(args.active_tol)
    smooth_weight = mp.mpf(args.smooth_weight)
    lo = mp.mpf(args.constraint_min)
    hi = mp.mpf(args.constraint_max)
    length = hi - lo
    max_q = args.jet_order - 1

    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(
        args,
        max(args.max_trace_q, args.basis - args.local_constraint_offset - 1),
    )
    base = local_case(args, args.basis, e_derivs_at_s0)
    source_nodes_used, source_rows = stack_source_rows(args, base["polys"], e_derivs_at_s0)
    source_on_scaled = source_rows * base["scaledModes"]
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], active_tol)
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis
    source_active = source_on_scaled * active_basis

    centers, R_global = trace_matrix_for_samples(args, base["polys"], args.trace_samples)
    weights = weights_for_centers(centers)
    trace_active = R_global * active_modes

    # Rank sanity check against the raw finite range model on the same grid.
    gram = trace_active.T * trace_active
    _gvals, _gkeep, gplus, _gtol = inverse_spd(gram, args.trace_tol)
    coeff_grid = source_active * gplus * trace_active.T
    grid_range_residual = source_active - coeff_grid * trace_active

    e_derivs_by_center, min_gap, min_cos, max_center_error = field_derivatives(
        args,
        centers,
        max_q,
    )
    kpolys = shifted_legendre_polys(args.kernel_degree, length)
    design = polynomial_design(kpolys, centers, weights, trace_active, lo)

    rows = []
    max_range_rel = mp.mpf("0")
    max_ibp_rel = mp.mpf("0")
    max_eta_l2 = mp.mpf("0")
    max_density_l2 = mp.mpf("0")
    for row_idx in range(source_active.rows):
        target = [source_active[row_idx, j] for j in range(source_active.cols)]
        coeffs = regularized_coefficients(
            design,
            target,
            args.kernel_degree,
            args.smooth_order,
            smooth_weight,
        )
        fit = [
            mp.fsum(design[j, m] * coeffs[m] for m in range(args.kernel_degree))
            for j in range(source_active.cols)
        ]
        range_defect = [fit[j] - target[j] for j in range(len(target))]
        range_rel = mp.sqrt(mp.fsum(x * x for x in range_defect)) / max(
            mp.mpf("1"),
            mp.sqrt(mp.fsum(x * x for x in target)),
        )

        kvals = kernel_values(kpolys, coeffs, centers, lo)
        eta_values = []
        for idx, center in enumerate(centers):
            k_derivs = kernel_derivatives(kpolys, coeffs, center, lo, max_q)
            eta_values.append(
                adjoint_source_value(e_derivs_by_center[idx], k_derivs, args.jet_order)
            )
        left_k = kernel_derivatives(kpolys, coeffs, centers[0], lo, max_q)
        right_k = kernel_derivatives(kpolys, coeffs, centers[-1], lo, max_q)
        left_b = trace_green_concomitant_row(e_derivs_by_center[0], left_k, args.jet_order)
        right_b = trace_green_concomitant_row(e_derivs_by_center[-1], right_k, args.jet_order)
        rhs_poly = integral_poly_row(base["polys"], centers, weights, eta_values)
        rhs_poly += boundary_poly_row(base["polys"], centers[-1], right_b)
        rhs_poly -= boundary_poly_row(base["polys"], centers[0], left_b)
        rhs_active = rhs_poly * active_modes
        lhs_active = mp.matrix(1, len(target))
        for j, value in enumerate(fit):
            lhs_active[0, j] = value
        ibp_defect = rhs_active - lhs_active
        ibp_rel = row_norm(ibp_defect) / max(mp.mpf("1"), row_norm(lhs_active))
        eta_l2 = mp.sqrt(mp.fsum(weights[i] * eta_values[i] ** 2 for i in range(len(weights))))
        density_l2 = mp.sqrt(mp.fsum(weights[i] * kvals[i] ** 2 for i in range(len(weights))))
        max_range_rel = max(max_range_rel, range_rel)
        max_ibp_rel = max(max_ibp_rel, ibp_rel)
        max_eta_l2 = max(max_eta_l2, eta_l2)
        max_density_l2 = max(max_density_l2, density_l2)
        component = "boundary" if row_idx % 2 == 0 else "adjointEval"
        rows.append(
            {
                "sourceNode": f(source_nodes_used[row_idx // 2]),
                "component": component,
                "rangeRelativeDefect": f(range_rel),
                "ibpRelativeDefectOnActive": f(ibp_rel),
                "etaL2Trap": f(eta_l2),
                "densityL2Trap": f(density_l2),
                "densityMaxAbs": f(max(abs(x) for x in kvals)),
                "coefficientL2": f(mp.sqrt(mp.fsum(c * c for c in coeffs))),
                "leftBoundaryRowNorm": f(list_norm(left_b)),
                "rightBoundaryRowNorm": f(list_norm(right_b)),
            }
        )

    grid_rel = frob_norm(grid_range_residual) / max(mp.mpf("1"), frob_norm(source_active))
    data = {
        "theoremName": "regularized finite adjoint Green BVP model",
        "basis": args.basis,
        "traceSamples": args.trace_samples,
        "kernelDegree": args.kernel_degree,
        "smoothOrder": args.smooth_order,
        "smoothWeight": f(smooth_weight),
        "interval": [f(lo), f(hi)],
        "activeDim": len(active_idx),
        "activeFloor": f(active_floor),
        "rawGridRangeResidualRelative": f(grid_rel),
        "regularizedRangeResidualRelative": f(max_range_rel),
        "regularizedMaxIbpRelativeDefectOnActive": f(max_ibp_rel),
        "regularizedMaxEtaL2": f(max_eta_l2),
        "regularizedMaxDensityL2": f(max_density_l2),
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "maxCenterMismatch": f(max_center_error),
        "rows": rows,
        "statuses": {
            "exactActiveRangeOnRegularizedKernel": status(
                "exact active range on regularized kernels",
                max_range_rel < mp.mpf("1e-40"),
                "Polynomial densities solve the two active moment equations.",
            ),
            "regularizedAdjointIbpSmall": status(
                "regularized adjoint IBP diagnostic",
                max_ibp_rel < mp.mpf("1e-3"),
                "Tests whether the smooth finite kernels behave like solutions of P^*K=eta with matched endpoint concomitant.",
            ),
            "continuumBvpClosed": status(
                "continuum adjoint BVP compactness",
                False,
                "Still requires a uniform norm bound and convergence of the regularized kernels in the trace-dual topology.",
            ),
        },
        "interpretation": (
            "This is the finite regularized replacement for differentiating the "
            "raw pseudoinverse trace coefficients.  A stable small IBP defect "
            "under degree/sample refinement would be the next numerical signal "
            "for the continuum adjoint Green proof."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Regularized adjoint Green BVP model", flush=True)
    print(f"  active dim: {len(active_idx)}", flush=True)
    print(f"  raw grid range residual: {fmt(grid_rel, 12)}", flush=True)
    print(f"  regularized range residual: {fmt(max_range_rel, 12)}", flush=True)
    print(f"  regularized max IBP defect: {fmt(max_ibp_rel, 12)}", flush=True)
    print(f"  regularized max density L2: {fmt(max_density_l2, 12)}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
