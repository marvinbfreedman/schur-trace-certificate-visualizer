#!/usr/bin/env python3
r"""Adjoint Green boundary diagnostic for the trace-to-source kernel.

The finite active range identity gives sampled coefficients C_N such that

    E_active = C_N R_N.

Interpreting C_{N,i}=w_i K_N(a_i), the Lagrange identity predicts

    int K_N Lambda_a(f) da
      = int f(a) P^*K_N(a) da + [B_P[K_N,f]]_{a_-}^{a_+}.

This script computes a finite-difference approximation of P^*K_N and the
endpoint concomitant, then checks this identity on the source-active high
block.  It is a diagnostic for the continuum boundary-value problem, not a
rigorous proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

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
from lambda_differential_closure import coefficient_derivatives, finite_diff_weights  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    poly_derivative_value,
    poly_value,
    trace_matrix,
)
from source_concomitant_membership import sample_trace_field  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def weights_for_centers(centers):
    if len(centers) == 1:
        return [mp.mpf("1")]
    weights = []
    for idx in range(len(centers)):
        if idx == 0:
            weights.append(abs(centers[1] - centers[0]) / 2)
        elif idx == len(centers) - 1:
            weights.append(abs(centers[-1] - centers[-2]) / 2)
        else:
            weights.append(abs(centers[idx + 1] - centers[idx - 1]) / 2)
    return weights


def active_indices(source_vals, source_top, active_tol):
    floor = active_tol * max(mp.mpf("1"), source_top)
    return [idx for idx, val in enumerate(source_vals) if val > floor], floor


def trace_args(args, constraints: int) -> SimpleNamespace:
    out = child_args(args, args.basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


def row_from_values(values):
    out = mp.matrix(1, len(values))
    for idx, value in enumerate(values):
        out[0, idx] = value
    return out


def density_derivatives(centers, density, idx: int, max_q: int, width: int):
    out = []
    for deriv in range(max_q + 1):
        stencil, fd_weights = finite_diff_weights(centers, idx, deriv, width)
        out.append(mp.fsum(fd_weights[j] * density[stencil[j]] for j in range(len(stencil))))
    return out


def boundary_poly_row(polys, point, brow):
    out = mp.matrix(1, len(polys))
    for col, poly in enumerate(polys):
        out[0, col] = mp.fsum(
            brow[j] * poly_derivative_value(poly, point, j)
            for j in range(len(brow))
        )
    return out


def integral_poly_row(polys, centers, weights, eta_values):
    out = mp.matrix(1, len(polys))
    for col, poly in enumerate(polys):
        out[0, col] = mp.fsum(
            weights[i] * eta_values[i] * poly_value(poly, centers[i])
            for i in range(len(centers))
        )
    return out


def row_norm(row):
    return mp.sqrt(mp.fsum(abs(row[0, j]) ** 2 for j in range(row.cols)))


def max_abs(values):
    return max((abs(value) for value in values), default=mp.mpf("0"))


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
    parser.add_argument("--stencil-width", type=int, default=9)
    parser.add_argument("--json-out", default="adjoint_green_boundary_diagnostic.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    active_tol = mp.mpf(args.active_tol)
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

    targs = trace_args(args, args.global_constraints)
    centers, R_global = trace_matrix(base["polys"], targs)
    weights = weights_for_centers(centers)
    trace_active = R_global * active_modes
    gram = trace_active.T * trace_active
    _gvals, _gkeep, gplus, _gtol = inverse_spd(gram, args.trace_tol)
    coeff = source_active * gplus * trace_active.T
    active_from_trace = coeff * trace_active
    active_residual = source_active - active_from_trace

    field_args = SimpleNamespace(**vars(args))
    field_args.s_min = args.constraint_min
    field_args.s_max = args.constraint_max
    field_args.samples = args.global_constraints
    centers_field, vecs, min_gap, min_cos = sample_trace_field(field_args)
    max_center_error = max(abs(centers[i] - centers_field[i]) for i in range(len(centers)))
    max_q = args.jet_order - 1
    width = min(args.stencil_width, len(centers))
    if width <= max_q:
        raise SystemExit("--stencil-width/global-constraints must exceed jet_order-1")

    e_derivs_by_center = [
        coefficient_derivatives(centers, vecs, idx, max_q, width)
        for idx in range(len(centers))
    ]

    rows = []
    max_ibp_rel = mp.mpf("0")
    max_active_rel = mp.mpf("0")
    max_eta_l2 = mp.mpf("0")
    for row_idx in range(coeff.rows):
        coeff_values = [coeff[row_idx, j] for j in range(coeff.cols)]
        density = [coeff_values[i] / weights[i] for i in range(len(weights))]
        eta_values = []
        for idx in range(len(centers)):
            k_derivs = density_derivatives(centers, density, idx, max_q, width)
            eta_values.append(
                adjoint_source_value(e_derivs_by_center[idx], k_derivs, args.jet_order)
            )

        left_k = density_derivatives(centers, density, 0, max_q, width)
        right_k = density_derivatives(centers, density, len(centers) - 1, max_q, width)
        left_b = trace_green_concomitant_row(e_derivs_by_center[0], left_k, args.jet_order)
        right_b = trace_green_concomitant_row(
            e_derivs_by_center[-1],
            right_k,
            args.jet_order,
        )
        lhs_poly = row_from_values(coeff_values) * R_global
        rhs_poly = integral_poly_row(base["polys"], centers, weights, eta_values)
        rhs_poly += boundary_poly_row(base["polys"], centers[-1], right_b)
        rhs_poly -= boundary_poly_row(base["polys"], centers[0], left_b)

        lhs_active = lhs_poly * active_modes
        rhs_active = rhs_poly * active_modes
        ibp_defect = rhs_active - lhs_active
        source_target = source_active[row_idx, :]
        active_defect = lhs_active - source_target
        lhs_norm = row_norm(lhs_active)
        source_norm = row_norm(source_target)
        ibp_rel = row_norm(ibp_defect) / max(mp.mpf("1"), lhs_norm)
        active_rel = row_norm(active_defect) / max(mp.mpf("1"), source_norm)
        eta_l2 = mp.sqrt(mp.fsum(weights[i] * eta_values[i] ** 2 for i in range(len(weights))))
        max_ibp_rel = max(max_ibp_rel, ibp_rel)
        max_active_rel = max(max_active_rel, active_rel)
        max_eta_l2 = max(max_eta_l2, eta_l2)
        component = "boundary" if row_idx % 2 == 0 else "adjointEval"
        rows.append(
            {
                "sourceNode": f(source_nodes_used[row_idx // 2]),
                "component": component,
                "ibpRelativeDefectOnActive": f(ibp_rel),
                "activeRangeRelativeDefect": f(active_rel),
                "etaL2Trap": f(eta_l2),
                "etaMaxAbs": f(max_abs(eta_values)),
                "densityL2Trap": f(
                    mp.sqrt(mp.fsum(weights[i] * density[i] ** 2 for i in range(len(weights))))
                ),
                "densityMaxAbs": f(max_abs(density)),
                "leftBoundaryRowNorm": f(mp.sqrt(mp.fsum(x * x for x in left_b))),
                "rightBoundaryRowNorm": f(mp.sqrt(mp.fsum(x * x for x in right_b))),
            }
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "globalConstraints": args.global_constraints,
        "interval": [f(mp.mpf(args.constraint_min)), f(mp.mpf(args.constraint_max))],
        "activeTol": f(active_tol),
        "activeDim": len(active_idx),
        "activeFloor": f(active_floor),
        "sourceTop": f(base["sourceTop"]),
        "traceCenters": [f(x) for x in centers],
        "traceWeights": [f(x) for x in weights],
        "maxCenterMismatch": f(max_center_error),
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "activeRangeFrobeniusRelative": f(
            frob_norm(active_residual) / max(mp.mpf("1"), frob_norm(source_active))
        ),
        "maxIbpRelativeDefectOnActive": f(max_ibp_rel),
        "maxActiveRangeRelativeDefect": f(max_active_rel),
        "maxEtaL2Trap": f(max_eta_l2),
        "rows": rows,
        "status": {
            "finiteActiveRangeClosed": max_active_rel < mp.mpf("1e-40"),
            "discreteAdjointIbpSmall": max_ibp_rel < mp.mpf("1e-3"),
            "continuumAdjointGreenClosed": False,
        },
        "interpretation": (
            "Finite-difference P^*K_N and endpoint concomitant diagnostic. "
            "Large IBP defects indicate the need for an analytic adjoint "
            "Green construction rather than differentiating coarse sampled "
            "kernel coefficients."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Adjoint Green boundary diagnostic")
    print(f"  active dim: {len(active_idx)}")
    print(f"  active range relative: {fmt(mp.mpf(data['activeRangeFrobeniusRelative']), 12)}")
    print(f"  max IBP relative defect: {fmt(max_ibp_rel, 12)}")
    print(f"  max eta L2: {fmt(max_eta_l2, 12)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
