#!/usr/bin/env python3
r"""Endpoint selection using exact local trace derivatives on each ODE segment.

``adjoint_green_endpoint_selection.py`` builds the adjoint companion system by
sampling the moving endpoint-defect eigenrow and finite-differencing those
samples.  Mesh-refinement tests showed that this is too unstable to serve as a
center for a rank-ball proof.

This script keeps the same endpoint map definition but replaces the sampled
finite-difference coefficient field by analytic confluent eigen-derivatives at
each ODE mesh point.  At a mesh point ``a`` it computes derivatives of the
trace row ``e(a)`` from the confluent eigen-equation for the endpoint kernel,
then forms the adjoint companion matrix for ``P^*K=0``.

The result is still not a rigorous interval enclosure.  It is a better center
candidate for the next validated-ODE step.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    endpoint_active_row,
    endpoint_map,
    left_obstruction,
    mat_frobenius,
    matrix_rank_gram,
    min_norm_solve,
    propagate_flows,
    row_norm,
    status,
)
from adjoint_green_jump_conditions import jump_matrix, solve_jump  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def mesh_points(left: mp.mpf, right: mp.mpf, samples: int) -> list[mp.mpf]:
    if samples < 2:
        raise SystemExit("--ode-samples must be at least 2")
    h = (right - left) / (samples - 1)
    return [left + i * h for i in range(samples)]


def exact_segment_derivatives(args, centers: list[mp.mpf], max_q: int):
    out = []
    rows = []
    previous = None
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    for idx, center in enumerate(centers):
        local_args = SimpleNamespace(**vars(args))
        local_args.s0 = str(center)
        vals, e_derivs, lam_derivs = exact_trace_derivatives(local_args, max_q)
        if previous is not None:
            cos = mp.fsum(previous[k] * e_derivs[0][k] for k in range(len(previous)))
            if cos < 0:
                e_derivs = [[-value for value in row] for row in e_derivs]
                cos = -cos
            min_cos = min(min_cos, abs(cos))
        previous = list(e_derivs[0])
        gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
        min_gap = min(min_gap, gap)
        out.append(e_derivs)
        rows.append(
            {
                "index": idx,
                "s": f(center),
                "lambda0": f(vals[0]),
                "lambda1": f(vals[1]),
                "gap": f(gap),
                "e8": f(e_derivs[0][-1]),
            }
        )
    return out, rows, min_gap, min_cos


def endpoint_rows_for_sources(args, base, active_modes, e_left, e_right, flow_right, e_s0):
    s0 = mp.mpf(args.s0)
    left = mp.mpf(args.constraint_min)
    right = mp.mpf(args.constraint_max)
    jmat = jump_matrix(e_s0, args.jet_order - 1)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    zero_left = mp.matrix(args.jet_order - 1, 1)
    rows = []
    for t in t_values:
        h_derivs = source_derivatives(mp.pi, s0, t, args.jet_order, r_nodes, r_weights)
        boundary_source = trace_green_concomitant_row(e_s0, h_derivs, args.jet_order)
        pstar = adjoint_source_value(e_s0, h_derivs, args.jet_order)
        eval_source = [pstar] + [mp.mpf("0") for _ in range(args.jet_order - 2)]
        for component, source_row in (("boundary", boundary_source), ("adjointEval", eval_source)):
            delta, jump_residual, _target = solve_jump(jmat, source_row)
            delta_vec = mp.matrix(len(delta), 1)
            for i, value in enumerate(delta):
                delta_vec[i] = value
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
            rows.append(
                {
                    "sourceNode": f(t),
                    "component": component,
                    "jumpNorm": f(row_norm(delta)),
                    "endpointVector": [f(value) for value in b_values],
                    "endpointBeforeNorm": f(row_norm(b_values)),
                    "jumpResidualNorm": f(row_norm(jump_residual)),
                }
            )
    return rows


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
    parser.add_argument("--ode-samples", type=int, default=7)
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
    parser.add_argument("--json-out", default="adjoint_green_endpoint_selection_exact_segments.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
    left = mp.mpf(args.constraint_min)
    right = mp.mpf(args.constraint_max)
    centers = mesh_points(left, right, args.ode_samples)
    idx0 = min(range(len(centers)), key=lambda idx: abs(centers[idx] - s0))
    if abs(centers[idx0] - s0) > mp.mpf("1e-30"):
        raise SystemExit("s0 must lie on the ODE sample grid")

    print(
        f"Exact-segment endpoint selection samples={args.ode_samples} "
        f"max_q={max(args.max_trace_q, args.jet_order - 1)}",
        flush=True,
    )
    derivs, derivative_rows, min_gap, min_cos = exact_segment_derivatives(
        args,
        centers,
        max(args.max_trace_q, args.jet_order - 1),
    )
    e_s0 = derivs[idx0]
    base, active_modes, active_idx = active_setup(args, e_s0)
    flow_left, flow_right, _idx0 = propagate_flows(centers, derivs, s0)
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
    rank, eigs, rank_threshold = matrix_rank_gram(
        endpoint_mat * endpoint_mat.T,
        mp.mpf(args.trace_tol),
    )
    obstruction, obstruction_eigs = left_obstruction(endpoint_mat)
    obstruction_residual = endpoint_mat.T * obstruction
    obstruction_residual_norm = row_norm(
        [obstruction_residual[i] for i in range(obstruction_residual.rows)]
    )
    source_rows = endpoint_rows_for_sources(
        args,
        base,
        active_modes,
        e_left,
        e_right,
        flow_right,
        e_s0,
    )
    max_residual = mp.mpf("0")
    max_z_norm = mp.mpf("0")
    max_obstruction_pairing = mp.mpf("0")
    max_obstruction_relative_pairing = mp.mpf("0")
    rows = []
    for row in source_rows:
        b_values = [mp.mpf(str(value)) for value in row["endpointVector"]]
        z, endpoint_residual, _gram = min_norm_solve(
            endpoint_mat,
            [-value for value in b_values],
        )
        residual_norm = row_norm([endpoint_residual[i] for i in range(endpoint_residual.rows)])
        z_norm = row_norm([z[i] for i in range(z.rows)])
        beta_before = row_norm(b_values)
        obstruction_pairing = mp.fsum(
            obstruction[i] * b_values[i] for i in range(len(b_values))
        )
        obstruction_relative = abs(obstruction_pairing) / max(mp.mpf("1"), beta_before)
        max_residual = max(max_residual, residual_norm)
        max_z_norm = max(max_z_norm, z_norm)
        max_obstruction_pairing = max(max_obstruction_pairing, abs(obstruction_pairing))
        max_obstruction_relative_pairing = max(
            max_obstruction_relative_pairing,
            obstruction_relative,
        )
        row.update(
            {
                "leftObstructionPairing": f(obstruction_pairing),
                "leftObstructionRelativePairing": f(obstruction_relative),
                "selectedZNorm": f(z_norm),
                "endpointResidualNorm": f(residual_norm),
                "endpointRelativeResidual": f(residual_norm / max(mp.mpf("1"), beta_before)),
                "selectedZ": [f(z[i]) for i in range(z.rows)],
            }
        )
        rows.append(row)

    data = {
        "theoremName": "adjoint Green endpoint active selection with exact segment derivatives",
        "basis": args.basis,
        "interval": [f(left), f(right)],
        "s0": f(s0),
        "odeSamples": args.ode_samples,
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
        "derivativeRows": derivative_rows,
        "rows": rows,
        "statuses": {
            "activeEndpointMapFullRowRank": status(
                "active endpoint map full row rank",
                rank == len(active_idx),
                "Exact segment-derivative center has full sampled active row rank.",
            ),
            "actualEndpointRhsInRange": status(
                "actual endpoint RHS in Range(M)",
                max_residual < mp.mpf("1e-20"),
                "For the actual source rows, the endpoint vectors b(d) are in the computed endpoint-map range.",
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
            "Finite endpoint diagnostic using analytic confluent eigen-derivatives "
            "at each segment point instead of finite differences of the moving "
            "eigenrow."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Exact-segment adjoint Green endpoint active selection", flush=True)
    print(f"  active dim: {len(active_idx)}", flush=True)
    print(f"  endpoint map rank: {rank}", flush=True)
    print(f"  endpoint map Frobenius: {fmt(mat_frobenius(endpoint_mat), 12)}", flush=True)
    print(f"  flow left/right Frobenius: {fmt(mat_frobenius(flow_left), 8)} / {fmt(mat_frobenius(flow_right), 8)}", flush=True)
    print(f"  max endpoint residual: {fmt(max_residual, 12)}", flush=True)
    print(f"  max selected z norm: {fmt(max_z_norm, 12)}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
