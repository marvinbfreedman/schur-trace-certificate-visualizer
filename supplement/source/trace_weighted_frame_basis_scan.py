#!/usr/bin/env python3
r"""Scan the weighted trace frame bound across Galerkin bases.

For the active source subspace, let

    T_N = W_N^{1/2} R_N,
    E_N = source-active rows.

The minimal discrete L2 density solving

    E_N = Y_N T_N

is ``Y_N = E_N T_N^+``.  Therefore

    ||Y_N|| <= ||E_N|| / s_min(T_N).

This is the finite frame inequality behind the continuum range theorem.  A
stable positive lower bound for ``T_N^*T_N`` on the active source subspace is
the numerical version of interval trace observability.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_active_range_inclusion import inverse_spd, op_norm  # noqa: E402
from global_trace_observability_gap import (  # noqa: E402
    child_args,
    f,
    fmt,
    frob_norm,
    parse_ints,
    stack_source_rows,
)
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402
from trace_to_source_kernel_refinement import (  # noqa: E402
    active_indices,
    scale_rows,
    weighted_density_stats,
    weights_for_centers,
)


def trace_args(args, basis: int, constraints: int):
    out = child_args(args, basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


def max_weighted_row_l2(weighted_y, weights, nodes):
    max_l2 = mp.mpf("0")
    max_l1 = mp.mpf("0")
    max_point = mp.mpf("0")
    rows = []
    for source_i, node in enumerate(nodes):
        for component in range(2):
            row = 2 * source_i + component
            stats = weighted_density_stats(
                [weighted_y[row, j] for j in range(weighted_y.cols)],
                weights,
            )
            max_l2 = max(max_l2, mp.mpf(stats["densityL2"]))
            max_l1 = max(max_l1, mp.mpf(stats["densityL1"]))
            max_point = max(max_point, mp.mpf(stats["maxDensity"]))
            rows.append(
                {
                    "sourceNode": f(node),
                    "component": "boundary" if component == 0 else "adjointEval",
                    **stats,
                }
            )
    return max_l2, max_l1, max_point, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20")
    parser.add_argument("--trace-counts", default="9 13")
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
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
    parser.add_argument("--json-out", default="trace_weighted_frame_basis_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    trace_counts = parse_ints(args.trace_counts)
    active_tol = mp.mpf(args.active_tol)
    max_local_constraints = max(1, max(bases) - args.local_constraint_offset)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    rows = []
    print(
        f"Weighted trace frame basis scan bases={bases} trace_counts={trace_counts}",
        flush=True,
    )
    print(
        "  basis traces active rank frame_min frame_max density_op max_row_l2 bound_op",
        flush=True,
    )

    for basis in bases:
        base = local_case(args, basis, e_derivs)
        nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs)
        source_on_scaled = source_rows * base["scaledModes"]
        active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], active_tol)
        active_basis = columns(base["sourceVecs"], active_idx)
        source_active = source_on_scaled * active_basis
        source_op = op_norm(source_active)
        source_frob = frob_norm(source_active)

        for count in trace_counts:
            centers, R_global = trace_matrix(base["polys"], trace_args(args, basis, count))
            weights = weights_for_centers(centers)
            sqrt_weights = [mp.sqrt(w) for w in weights]
            trace_active = R_global * base["scaledModes"] * active_basis
            weighted_trace = scale_rows(trace_active, sqrt_weights)
            frame = weighted_trace.T * weighted_trace
            frame_vals, frame_keep, frame_plus, frame_tol = inverse_spd(frame, args.trace_tol)
            weighted_y = source_active * frame_plus * weighted_trace.T
            residual = source_active - weighted_y * weighted_trace
            residual_norm = frob_norm(residual)

            frame_min = min(frame_vals) if len(frame_vals) else mp.mpf("0")
            frame_max = max(frame_vals) if len(frame_vals) else mp.mpf("0")
            density_op = op_norm(weighted_y)
            density_frob = frob_norm(weighted_y)
            max_row_l2, max_row_l1, max_point, row_stats = max_weighted_row_l2(
                weighted_y,
                weights,
                nodes,
            )
            bound_op = source_op / mp.sqrt(frame_min) if frame_min > 0 else mp.inf
            rows.append(
                {
                    "basis": basis,
                    "traceCount": count,
                    "traceMin": f(centers[0]) if centers else None,
                    "traceMax": f(centers[-1]) if centers else None,
                    "activeTol": f(active_tol),
                    "activeFloor": f(active_floor),
                    "activeDim": len(active_idx),
                    "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
                    "sourceTop": f(base["sourceTop"]),
                    "sourceOperator": f(source_op),
                    "sourceFrobenius": f(source_frob),
                    "frameRank": len(frame_keep),
                    "frameEigenvalues": [f(x) for x in frame_vals],
                    "frameTol": f(frame_tol),
                    "frameMin": f(frame_min),
                    "frameMax": f(frame_max),
                    "densityOperator": f(density_op),
                    "densityFrobenius": f(density_frob),
                    "maxRowDensityL2": f(max_row_l2),
                    "maxRowDensityL1": f(max_row_l1),
                    "maxPointDensity": f(max_point),
                    "frameBoundOperator": f(bound_op),
                    "boundRatio": f(density_op / bound_op if bound_op else mp.mpf("0")),
                    "rangeResidualRelative": f(residual_norm / source_frob if source_frob else mp.mpf("0")),
                    "weightedRowStats": row_stats,
                }
            )
            print(
                f"  {basis:5d} {count:6d} {len(active_idx):6d} {len(frame_keep):4d} "
                f"{fmt(frame_min, 8):>9} {fmt(frame_max, 8):>9} "
                f"{fmt(density_op, 8):>10} {fmt(max_row_l2, 8):>10} "
                f"{fmt(bound_op, 8):>10}",
                flush=True,
            )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceGrid": args.source_grid,
        "bases": bases,
        "traceCounts": trace_counts,
        "constraintMin": f(mp.mpf(args.constraint_min)),
        "constraintMax": f(mp.mpf(args.constraint_max)),
        "cutoff": args.cutoff,
        "activeTol": f(active_tol),
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Finite weighted frame inequality on the source-active subspace. "
            "Uniform lower bounds for frameMin imply uniform L2 bounds for "
            "the minimal Green density via ||Y|| <= ||E||/sqrt(frameMin)."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
