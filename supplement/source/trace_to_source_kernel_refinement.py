#!/usr/bin/env python3
r"""Refine the sampled trace-to-source kernel and test continuum scaling.

If

    E_active f = sum_j C_j Lambda_{a_j}(f)

is the quadrature approximation of

    E_active f = int_I K(a) Lambda_a(f) da,

then C_j should behave like w_j K(a_j).  This script scans the number of trace
samples and reports both the coefficient norms and the quadrature-scaled density
norms

    ||K||_L1 ~ sum_j |C_j|,
    ||K||_L2 ~ (sum_j |C_j|^2 / w_j)^(1/2).

Bounded density norms support an L1/L2 Green family.  Growth indicates that the
continuum range element may live in a rougher trace-dual space.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_active_range_inclusion import inverse_spd, op_norm  # noqa: E402
from global_trace_observability_gap import (  # noqa: E402
    child_args,
    f,
    fmt,
    frob_norm,
    source_nodes,
    stack_source_rows,
)
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402


def weights_for_centers(centers):
    if len(centers) == 1:
        return [mp.mpf("1")]
    weights = []
    for i in range(len(centers)):
        if i == 0:
            weights.append(abs(centers[1] - centers[0]) / 2)
        elif i == len(centers) - 1:
            weights.append(abs(centers[-1] - centers[-2]) / 2)
        else:
            weights.append(abs(centers[i + 1] - centers[i - 1]) / 2)
    return weights


def row_density_stats(coeffs, weights):
    vals = [mp.mpf(v) for v in coeffs]
    l1_density = mp.fsum(abs(v) for v in vals)
    l2_density = mp.sqrt(mp.fsum(abs(vals[i]) ** 2 / weights[i] for i in range(len(vals))))
    max_density = max(abs(vals[i]) / weights[i] for i in range(len(vals))) if vals else mp.mpf("0")
    max_coeff = max(abs(v) for v in vals) if vals else mp.mpf("0")
    variation_coeff = mp.fsum(abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1))
    return {
        "coefficientL1": f(l1_density),
        "densityL2": f(l2_density),
        "maxCoefficient": f(max_coeff),
        "maxDensity": f(max_density),
        "coefficientVariation": f(variation_coeff),
    }


def weighted_density_stats(yvals, weights):
    """Stats for y_j=sqrt(w_j) K(a_j), the discrete L2 density coordinates."""
    vals = [mp.mpf(v) for v in yvals]
    density_l2 = mp.sqrt(mp.fsum(abs(v) ** 2 for v in vals))
    density_l1 = mp.fsum(mp.sqrt(weights[i]) * abs(vals[i]) for i in range(len(vals)))
    max_density = max(abs(vals[i]) / mp.sqrt(weights[i]) for i in range(len(vals))) if vals else mp.mpf("0")
    max_coefficient = max(mp.sqrt(weights[i]) * abs(vals[i]) for i in range(len(vals))) if vals else mp.mpf("0")
    variation_density_coord = mp.fsum(abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1))
    return {
        "densityL1": f(density_l1),
        "densityL2": f(density_l2),
        "maxDensity": f(max_density),
        "maxCoefficient": f(max_coefficient),
        "densityCoordinateVariation": f(variation_density_coord),
    }


def scale_rows(mat, scales):
    out = mp.matrix(mat.rows, mat.cols)
    for i, scale in enumerate(scales):
        for j in range(mat.cols):
            out[i, j] = scale * mat[i, j]
    return out


def active_indices(source_vals, source_top, active_tol):
    floor = active_tol * max(mp.mpf("1"), source_top)
    return [i for i, val in enumerate(source_vals) if val > floor], floor


def trace_args(args, constraints: int) -> SimpleNamespace:
    out = child_args(args, args.basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


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
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--trace-counts", default="7 9 11 13")
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
    parser.add_argument("--json-out", default="trace_to_source_kernel_refinement.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    active_tol = mp.mpf(args.active_tol)
    trace_counts = [int(piece) for piece in args.trace_counts.replace(",", " ").split()]
    max_local_constraints = max(1, args.basis - args.local_constraint_offset)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    base = local_case(args, args.basis, e_derivs)
    nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs)
    source_on_scaled = source_rows * base["scaledModes"]
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], active_tol)
    active_basis = columns(base["sourceVecs"], active_idx)
    source_active = source_on_scaled * active_basis

    rows = []
    print(
        f"Trace-to-source kernel refinement basis={args.basis} active_dim={len(active_idx)}",
        flush=True,
    )
    print(
        "  traces rank rel_res coeff_op raw_L2 weighted_L2 weighted_L1 weighted_max",
        flush=True,
    )

    for count in trace_counts:
        centers, R_global = trace_matrix(base["polys"], trace_args(args, count))
        weights = weights_for_centers(centers)
        trace_active = R_global * base["scaledModes"] * active_basis
        gram = trace_active.T * trace_active
        gvals, gkeep, gplus, gtol = inverse_spd(gram, args.trace_tol)
        coeff = source_active * gplus * trace_active.T
        residual = source_active - coeff * trace_active
        source_norm = frob_norm(source_active)
        residual_norm = frob_norm(residual)

        sqrt_weights = [mp.sqrt(w) for w in weights]
        weighted_trace = scale_rows(trace_active, sqrt_weights)
        weighted_gram = weighted_trace.T * weighted_trace
        wgvals, wgkeep, wgplus, wgtol = inverse_spd(weighted_gram, args.trace_tol)
        weighted_y = source_active * wgplus * weighted_trace.T
        weighted_residual = source_active - weighted_y * weighted_trace
        weighted_residual_norm = frob_norm(weighted_residual)

        row_stats = []
        weighted_row_stats = []
        for source_i, node in enumerate(nodes):
            for component in range(2):
                row = 2 * source_i + component
                coeffs = [coeff[row, j] for j in range(coeff.cols)]
                yvals = [weighted_y[row, j] for j in range(weighted_y.cols)]
                row_stats.append(
                    {
                        "sourceNode": f(node),
                        "component": "boundary" if component == 0 else "adjointEval",
                        **row_density_stats(coeffs, weights),
                    }
                )
                weighted_row_stats.append(
                    {
                        "sourceNode": f(node),
                        "component": "boundary" if component == 0 else "adjointEval",
                        **weighted_density_stats(yvals, weights),
                    }
                )

        max_l1 = max(mp.mpf(row["coefficientL1"]) for row in row_stats) if row_stats else mp.mpf("0")
        max_l2 = max(mp.mpf(row["densityL2"]) for row in row_stats) if row_stats else mp.mpf("0")
        max_density = max(mp.mpf(row["maxDensity"]) for row in row_stats) if row_stats else mp.mpf("0")
        weighted_max_l1 = max(mp.mpf(row["densityL1"]) for row in weighted_row_stats) if weighted_row_stats else mp.mpf("0")
        weighted_max_l2 = max(mp.mpf(row["densityL2"]) for row in weighted_row_stats) if weighted_row_stats else mp.mpf("0")
        weighted_max_density = max(mp.mpf(row["maxDensity"]) for row in weighted_row_stats) if weighted_row_stats else mp.mpf("0")
        row = {
            "basis": args.basis,
            "traceCount": count,
            "traceMin": f(centers[0]) if centers else None,
            "traceMax": f(centers[-1]) if centers else None,
            "traceStep": f(centers[1] - centers[0]) if len(centers) > 1 else None,
            "activeDim": len(active_idx),
            "traceRankOnActive": len(gkeep),
            "traceGramEigenvalues": [f(x) for x in gvals],
            "traceGramTol": f(gtol),
            "coefficientOperator": f(op_norm(coeff)),
            "coefficientFrobenius": f(frob_norm(coeff)),
            "sourceActiveFrobenius": f(source_norm),
            "rangeResidualRelative": f(residual_norm / source_norm if source_norm else mp.mpf("0")),
            "maxCoefficientL1AsDensityL1": f(max_l1),
            "maxDensityL2": f(max_l2),
            "maxPointDensity": f(max_density),
            "rowStats": row_stats,
            "weightedTraceGramEigenvalues": [f(x) for x in wgvals],
            "weightedTraceGramTol": f(wgtol),
            "weightedTraceRankOnActive": len(wgkeep),
            "weightedDensityCoordinateOperator": f(op_norm(weighted_y)),
            "weightedDensityCoordinateFrobenius": f(frob_norm(weighted_y)),
            "weightedRangeResidualRelative": f(
                weighted_residual_norm / source_norm if source_norm else mp.mpf("0")
            ),
            "weightedMaxDensityL1": f(weighted_max_l1),
            "weightedMaxDensityL2": f(weighted_max_l2),
            "weightedMaxPointDensity": f(weighted_max_density),
            "weightedRowStats": weighted_row_stats,
        }
        rows.append(row)
        print(
            f"  {count:6d} {len(gkeep):4d} "
            f"{fmt(residual_norm / source_norm if source_norm else 0, 8):>10} "
            f"{fmt(op_norm(coeff), 8):>10} {fmt(max_l2, 8):>7} "
            f"{fmt(weighted_max_l2, 8):>11} {fmt(weighted_max_l1, 8):>11} "
            f"{fmt(weighted_max_density, 8):>12}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "cutoff": args.cutoff,
        "localConstraints": base["localConstraints"],
        "sourceNodes": [f(x) for x in nodes],
        "constraintMin": f(mp.mpf(args.constraint_min)),
        "constraintMax": f(mp.mpf(args.constraint_max)),
        "activeTol": f(active_tol),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "sourceTop": f(base["sourceTop"]),
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Quadrature scaling scan for the trace-to-source kernel.  "
            "The unweighted coefficient vector C_j represents weights times "
            "one possible continuum density K(a_j).  The weighted fields use "
            "the correct minimal discrete L2 density coordinates y_j=sqrt(w_j)K(a_j), "
            "computed from the weighted trace Gram R^* W R."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
