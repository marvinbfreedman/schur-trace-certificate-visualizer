#!/usr/bin/env python3
r"""Endpoint selection for the adjoint Green free vector z.

After the interior source jump is fixed, every adjoint Green coefficient has
the form

    Y(a)=Phi(a,s0) z,                 a < s0,
    Y(a)=Phi(a,s0)(z+Delta(d)),       a > s0,

where Y=(K,K',...,K^(7)) and Delta(d) is the jump vector computed by
``adjoint_green_jump_conditions.py``.  The remaining endpoint concomitant on
the active source space is affine in z:

    beta_active(z;d) = M z + b(d).

This script constructs a finite diagnostic for that endpoint map and chooses
    the minimum-norm solution

    z = - M^T (M M^T)^(-1) b(d)

when the endpoint right-hand side lies in Range(M).  Full row rank of M is
sufficient but not necessary.  The homogeneous flow Phi is approximated from
the sampled moving trace field; the algebraic range condition itself is the
continuum theorem target.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_jump_conditions import jump_matrix, solve_jump  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, stack_source_rows  # noqa: E402
from lambda_differential_closure import coefficient_derivatives  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, poly_derivative_value  # noqa: E402
from source_concomitant_membership import sample_trace_field, source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def row_norm(values) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in values))


def mat_frobenius(mat: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def coefficient_derivative(e_derivs, deriv: int, k: int) -> mp.mpf:
    return e_derivs[deriv][k] / mp.factorial(k)


def adjoint_coefficients(e_derivs, order: int = 8):
    coeffs = []
    for m in range(order + 1):
        total = mp.mpf("0")
        for k in range(m, order + 1):
            total += (
                ((-1) ** k)
                * mp.binomial(k, m)
                * coefficient_derivative(e_derivs, k - m, k)
            )
        coeffs.append(total)
    return coeffs


def companion_matrix(e_derivs, order: int = 8) -> mp.matrix:
    coeffs = adjoint_coefficients(e_derivs, order)
    lead = coeffs[order]
    mat = mp.matrix(order, order)
    for i in range(order - 1):
        mat[i, i + 1] = 1
    for j in range(order):
        mat[order - 1, j] = -coeffs[j] / lead
    return mat


def sample_derivatives(args, max_q: int):
    field_args = SimpleNamespace(**vars(args))
    field_args.s_min = args.constraint_min
    field_args.s_max = args.constraint_max
    field_args.samples = args.ode_samples
    centers, vecs, min_gap, min_cos = sample_trace_field(field_args)
    derivs = [
        coefficient_derivatives(centers, vecs, idx, max_q, args.stencil_width)
        for idx in range(len(centers))
    ]
    return centers, derivs, min_gap, min_cos


def nearest_index(centers, value):
    return min(range(len(centers)), key=lambda idx: abs(centers[idx] - value))


def propagate_flows(centers, derivs, s0):
    order = len(derivs[0][0]) - 1
    idx0 = nearest_index(centers, s0)
    if abs(centers[idx0] - s0) > mp.mpf("1e-30"):
        raise SystemExit("s0 must lie on the ODE sample grid for this diagnostic")

    ident = mp.eye(order)
    right = ident.copy()
    for idx in range(idx0, len(centers) - 1):
        h = centers[idx + 1] - centers[idx]
        amat = companion_matrix(derivs[idx], order)
        right = mp.expm(h * amat) * right

    left = ident.copy()
    for idx in range(idx0, 0, -1):
        h = centers[idx - 1] - centers[idx]
        amat = companion_matrix(derivs[idx], order)
        left = mp.expm(h * amat) * left
    return left, right, idx0


def boundary_row_on_polys(polys, point, brow):
    out = mp.matrix(1, len(polys))
    for col, poly in enumerate(polys):
        out[0, col] = mp.fsum(
            brow[j] * poly_derivative_value(poly, point, j)
            for j in range(len(brow))
        )
    return out


def endpoint_active_row(polys, active_modes, e_left, e_right, y_left, y_right, left, right):
    brow_left = trace_green_concomitant_row(e_left, [y_left[i] for i in range(y_left.rows)], 9)
    brow_right = trace_green_concomitant_row(e_right, [y_right[i] for i in range(y_right.rows)], 9)
    row = boundary_row_on_polys(polys, right, brow_right)
    row -= boundary_row_on_polys(polys, left, brow_left)
    active = row * active_modes
    return [active[0, j] for j in range(active.cols)]


def endpoint_map(polys, active_modes, e_left, e_right, flow_left, flow_right, left, right):
    order = flow_left.rows
    mat = mp.matrix(active_modes.cols, order)
    zero = mp.matrix(order, 1)
    for col in range(order):
        unit = mp.matrix(order, 1)
        unit[col] = 1
        values = endpoint_active_row(
            polys,
            active_modes,
            e_left,
            e_right,
            flow_left * unit,
            flow_right * unit,
            left,
            right,
        )
        for row, value in enumerate(values):
            mat[row, col] = value
    base_zero = endpoint_active_row(
        polys,
        active_modes,
        e_left,
        e_right,
        zero,
        zero,
        left,
        right,
    )
    return mat, base_zero


def min_norm_solve(full_row_rank: mp.matrix, rhs: list[mp.mpf]):
    rhs_vec = mp.matrix(len(rhs), 1)
    for i, value in enumerate(rhs):
        rhs_vec[i] = value
    gram = full_row_rank * full_row_rank.T
    alpha = mp.lu_solve(gram, rhs_vec)
    sol = full_row_rank.T * alpha
    residual = full_row_rank * sol - rhs_vec
    return sol, residual, gram


def left_obstruction(mat: mp.matrix) -> tuple[mp.matrix, list[mp.mpf]]:
    gram = (mat * mat.T + mat * mat.T) / 2
    vals, vecs = mp.eigsy(gram, eigvals_only=False)
    idx = min(range(len(vals)), key=lambda i: abs(vals[i]))
    out = mp.matrix(mat.rows, 1)
    for i in range(mat.rows):
        out[i] = vecs[i, idx]
    nrm = row_norm([out[i] for i in range(out.rows)])
    if nrm:
        out /= nrm
    if out.rows and out[0] < 0:
        out *= -1
    return out, [vals[i] for i in range(len(vals))]


def matrix_rank_gram(gram: mp.matrix, tol: mp.mpf):
    vals = mp.eigsy((gram + gram.T) / 2, eigvals_only=True)
    scale = max([abs(v) for v in vals] + [mp.mpf("1")])
    threshold = tol * scale
    return sum(1 for v in vals if v > threshold), vals, threshold


def active_setup(args, e_derivs_at_s0):
    base = local_case(args, args.basis, e_derivs_at_s0)
    active_tol = mp.mpf(args.active_tol)
    active_idx = [
        idx
        for idx, val in enumerate(base["sourceVals"])
        if val > active_tol * max(mp.mpf("1"), base["sourceTop"])
    ]
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis
    return base, active_modes, active_idx


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
    parser.add_argument("--ode-samples", type=int, default=11)
    parser.add_argument("--stencil-width", type=int, default=9)
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
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
    parser.add_argument("--json-out", default="adjoint_green_endpoint_selection.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
    left = mp.mpf(args.constraint_min)
    right = mp.mpf(args.constraint_max)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, args.max_trace_q)
    base, active_modes, active_idx = active_setup(args, e_derivs_at_s0)

    centers, derivs, min_gap, min_cos = sample_derivatives(args, args.jet_order - 1)
    flow_left, flow_right, idx0 = propagate_flows(centers, derivs, s0)
    e_left = derivs[0]
    e_right = derivs[-1]
    endpoint_mat, _zero = endpoint_map(
        base["polys"],
        active_modes,
        e_left,
        e_right,
        flow_left,
        flow_right,
        left,
        right,
    )
    rank, eigs, rank_threshold = matrix_rank_gram(endpoint_mat * endpoint_mat.T, mp.mpf(args.trace_tol))
    obstruction, obstruction_eigs = left_obstruction(endpoint_mat)
    obstruction_residual = endpoint_mat.T * obstruction
    obstruction_residual_norm = row_norm(
        [obstruction_residual[i] for i in range(obstruction_residual.rows)]
    )

    jmat = jump_matrix(e_derivs_at_s0, args.jet_order - 1)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    r_nodes, r_weights = __import__("quotient_factorization_mp").endpoint_b_quadrature(
        mp.mpf(args.kernel_rmax),
        args.kernel_order,
    )

    rows = []
    max_residual = mp.mpf("0")
    max_z_norm = mp.mpf("0")
    max_obstruction_pairing = mp.mpf("0")
    max_obstruction_relative_pairing = mp.mpf("0")
    for t in t_values:
        h_derivs = source_derivatives(mp.pi, s0, t, args.jet_order, r_nodes, r_weights)
        boundary_source = trace_green_concomitant_row(e_derivs_at_s0, h_derivs, args.jet_order)
        pstar = adjoint_source_value(e_derivs_at_s0, h_derivs, args.jet_order)
        eval_source = [pstar] + [mp.mpf("0") for _ in range(args.jet_order - 2)]
        for component, source_row in (("boundary", boundary_source), ("adjointEval", eval_source)):
            delta, jump_residual, _target = solve_jump(jmat, source_row)
            delta_vec = mp.matrix(len(delta), 1)
            for i, value in enumerate(delta):
                delta_vec[i] = value
            zero_left = mp.matrix(len(delta), 1)
            b_values = endpoint_active_row(
                base["polys"],
                active_modes,
                e_left,
                e_right,
                zero_left,
                flow_right * delta_vec,
                left,
                right,
            )
            z, endpoint_residual, _gram = min_norm_solve(
                endpoint_mat,
                [-value for value in b_values],
            )
            residual_norm = mp.sqrt(mp.fsum(endpoint_residual[i] ** 2 for i in range(endpoint_residual.rows)))
            z_norm = mp.sqrt(mp.fsum(z[i] ** 2 for i in range(z.rows)))
            beta_before = row_norm(b_values)
            obstruction_pairing = mp.fsum(
                obstruction[i] * b_values[i] for i in range(len(b_values))
            )
            obstruction_relative = abs(obstruction_pairing) / max(
                mp.mpf("1"), beta_before
            )
            max_residual = max(max_residual, residual_norm)
            max_z_norm = max(max_z_norm, z_norm)
            max_obstruction_pairing = max(
                max_obstruction_pairing, abs(obstruction_pairing)
            )
            max_obstruction_relative_pairing = max(
                max_obstruction_relative_pairing, obstruction_relative
            )
            rows.append(
                {
                    "sourceNode": f(t),
                    "component": component,
                    "jumpNorm": f(row_norm(delta)),
                    "endpointVector": [f(value) for value in b_values],
                    "endpointBeforeNorm": f(beta_before),
                    "leftObstructionPairing": f(obstruction_pairing),
                    "leftObstructionRelativePairing": f(obstruction_relative),
                    "selectedZNorm": f(z_norm),
                    "endpointResidualNorm": f(residual_norm),
                    "endpointRelativeResidual": f(residual_norm / max(mp.mpf("1"), beta_before)),
                    "selectedZ": [f(z[i]) for i in range(z.rows)],
                }
            )

    data = {
        "theoremName": "adjoint Green endpoint active selection",
        "basis": args.basis,
        "interval": [f(left), f(right)],
        "s0": f(s0),
        "odeSamples": args.ode_samples,
        "stencilWidth": args.stencil_width,
        "activeDim": len(active_idx),
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "endpointMapRank": rank,
        "endpointMapGramEigenvalues": [f(value) for value in eigs],
        "endpointMapRankThreshold": f(rank_threshold),
        "endpointMap": [
            [f(endpoint_mat[i, j]) for j in range(endpoint_mat.cols)]
            for i in range(endpoint_mat.rows)
        ],
        "leftObstructionVector": [f(obstruction[i]) for i in range(obstruction.rows)],
        "leftObstructionGramEigenvalues": [f(value) for value in obstruction_eigs],
        "leftObstructionResidualNorm": f(obstruction_residual_norm),
        "maxLeftObstructionPairing": f(max_obstruction_pairing),
        "maxLeftObstructionRelativePairing": f(max_obstruction_relative_pairing),
        "endpointMapFrobenius": f(mat_frobenius(endpoint_mat)),
        "flowLeftFrobenius": f(mat_frobenius(flow_left)),
        "flowRightFrobenius": f(mat_frobenius(flow_right)),
        "maxEndpointResidualNorm": f(max_residual),
        "maxSelectedZNorm": f(max_z_norm),
        "rows": rows,
        "statuses": {
            "activeEndpointMapFullRowRank": status(
                "active endpoint map full row rank",
                rank == len(active_idx),
                "Sufficient but not necessary; the diagnostic may still close if the actual endpoint right-hand sides lie in Range(M).",
            ),
            "actualEndpointRhsInRange": status(
                "actual endpoint RHS in Range(M)",
                max_residual < mp.mpf("1e-20"),
                "For the actual source rows, the endpoint vectors b(d) are in the finite endpoint-map range.",
            ),
            "activeEndpointKilled": status(
                "active endpoint concomitant killed",
                max_residual < mp.mpf("1e-20"),
                "Minimum-norm z solves Mz=-b(d) on the active two-dimensional source space.",
            ),
            "finiteLeftObstructionAnnihilatesActualRows": status(
                "finite left obstruction annihilates actual rows",
                max_obstruction_relative_pairing < mp.mpf("1e-20"),
                "The resolved left endpoint obstruction annihilates the actual source endpoint vectors.",
            ),
            "uniformTraceDualBoundClosed": status(
                "uniform trace-dual bound",
                False,
                "Still needs a rigorous continuum bound for the exact fundamental matrix and M^+ over the source window.",
            ),
        },
        "interpretation": (
            "Finite diagnostic for the endpoint selection theorem.  The exact "
            "continuum proof needs the same range condition and uniform bounds "
            "for the true fundamental matrix/M^+, not the sampled ODE surrogate."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Adjoint Green endpoint active selection", flush=True)
    print(f"  active dim: {len(active_idx)}", flush=True)
    print(f"  endpoint map rank: {rank}", flush=True)
    print(f"  endpoint map Frobenius: {fmt(mat_frobenius(endpoint_mat), 12)}", flush=True)
    print(f"  flow left/right Frobenius: {fmt(mat_frobenius(flow_left), 8)} / {fmt(mat_frobenius(flow_right), 8)}", flush=True)
    print(f"  max endpoint residual: {fmt(max_residual, 12)}", flush=True)
    print(f"  max selected z norm: {fmt(max_z_norm, 12)}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
