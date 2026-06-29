#!/usr/bin/env python3
r"""Weighted source-quadrature refinement for source-side noncollapse.

The source-side rank theorem ultimately lives in a continuum source Hilbert
space over the source window.  The earlier source-side diagnostics used an
unweighted fixed source sample.  This script replaces that by trapezoid-weighted
source rows, so that the finite source matrix approximates

    E f(u) in L^2([u_min,u_max]; R^2).

It reports whether the active dimension, spectral gap, and rank margin survive
as the source grid is refined.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import eigvals_gram, parse_orders, submatrix_rows, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints, source_nodes, stack_source_rows  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from trace_active_derivative_rank import derivative_matrix, normalized_columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def weights_for_nodes(nodes):
    if len(nodes) == 1:
        return [mp.mpf("1")]
    weights = []
    for i in range(len(nodes)):
        if i == 0:
            weights.append(abs(nodes[1] - nodes[0]) / 2)
        elif i == len(nodes) - 1:
            weights.append(abs(nodes[-1] - nodes[-2]) / 2)
        else:
            weights.append(abs(nodes[i + 1] - nodes[i - 1]) / 2)
    return weights


def weight_source_rows(rows, weights):
    out = mp.matrix(rows.rows, rows.cols)
    for i in range(rows.rows):
        scale = mp.sqrt(weights[i // 2])
        for j in range(rows.cols):
            out[i, j] = scale * rows[i, j]
    return out


def columns(mat, indices):
    out = mp.matrix(mat.rows, len(indices))
    for j, idx in enumerate(indices):
        for i in range(mat.rows):
            out[i, j] = mat[i, idx]
    return out


def op_norm(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat * mat.T), eigvals_only=True)
    return mp.sqrt(max([mp.mpf("0"), *vals]))


def singular_values(mat):
    vals = eigvals_gram(mat)
    return [mp.sqrt(max(mp.mpf("0"), val)) for val in vals]


def grid_case(args, base, lmat, e_derivs_at_s0, grid):
    gargs = SimpleNamespace(**vars(args))
    gargs.source_grid = grid
    nodes, source_rows = stack_source_rows(gargs, base["polys"], e_derivs_at_s0)
    weights = weights_for_nodes(nodes)
    weighted_rows = weight_source_rows(source_rows, weights)
    source_on_scaled = weighted_rows * base["scaledModes"]
    source_side = sym(source_on_scaled * source_on_scaled.T)
    vals, vecs = mp.eigsy(source_side, eigvals_only=False)
    source_top = max([abs(v) for v in vals] + [mp.mpf("1")])
    active_idx, active_floor = active_indices(vals, source_top, mp.mpf(args.active_tol))
    row = {
        "basis": args.basis,
        "sourceGrid": grid,
        "sourceSideDim": source_side.rows,
        "sourceTop": f(source_top),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "sourceEigenvaluesTopDown": [f(vals[i]) for i in reversed(range(max(0, len(vals) - 6), len(vals)))],
    }
    if len(active_idx) != 2:
        row["interpretation"] = "active dimension is not two at this source quadrature"
        return row

    active_left = columns(vecs, active_idx)
    bmat = lmat * source_on_scaled.T
    compressed = bmat * active_left
    raw_svals = singular_values(compressed)
    normalized_response, response_norms = normalized_columns(compressed)
    normalized_svals = singular_values(normalized_response)
    b_norm = op_norm(bmat)

    active_min_idx = active_idx[0]
    next_below = vals[active_min_idx - 1] if active_min_idx > 0 else mp.mpf("0")
    active_min = vals[active_min_idx]
    spectral_gap = active_min - next_below
    raw_min = min(raw_svals)
    relative_rank_margin = raw_min / b_norm if b_norm else mp.mpf("0")
    gap_normed = spectral_gap / source_top if source_top else mp.mpf("0")
    row.update(
        {
            "activeSourceEigenvalues": [f(vals[i]) for i in active_idx],
            "activeSourceEigenvaluesFracOfTop": [f(vals[i] / source_top) for i in active_idx],
            "nextBelowActiveEigenvalue": f(next_below),
            "spectralGapToComplement": f(spectral_gap),
            "spectralGapToComplementFracOfTop": f(gap_normed),
            "rawCompressedSingularValues": [f(x) for x in raw_svals],
            "rawRankMargin": f(raw_min),
            "responseOperatorNorm": f(b_norm),
            "relativeRankMargin": f(relative_rank_margin),
            "responseColumnNorms": [f(x) for x in response_norms],
            "normalizedResponseSingularValues": [f(x) for x in normalized_svals],
            "normalizedRankMargin": f(min(normalized_svals)),
            "relativeBErrorMax": f(relative_rank_margin / 4),
            "topNormalizedTErrorMax": f(relative_rank_margin * gap_normed / 16),
        }
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--source-grids", default="5 7 9 11 13 17")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--max-deriv", type=int, default=8)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-6")
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
    parser.add_argument("--json-out", default="source_side_quadrature_refinement.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    grids = parse_ints(args.source_grids)
    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs_at_s0)

    point_args = SimpleNamespace(**vars(args))
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(point_args, base, base["scaledModes"], mp.mpf(args.point), max(args.orders))
    lmat = submatrix_rows(full_response, args.orders)

    print(
        f"Weighted source quadrature refinement basis={args.basis} point={args.point} orders={args.orders}",
        flush=True,
    )
    print("  grid dim eta rel_margin gap/top allow_B_rel allow_T/top", flush=True)
    rows = []
    for grid in grids:
        row = grid_case(args, base, lmat, e_derivs_at_s0, grid)
        rows.append(row)
        if row["activeDim"] == 2:
            print(
                f"  {grid:4d} {row['sourceSideDim']:3d} "
                f"{fmt(mp.mpf(row['normalizedRankMargin']), 8):>10} "
                f"{fmt(mp.mpf(row['relativeRankMargin']), 8):>10} "
                f"{fmt(mp.mpf(row['spectralGapToComplementFracOfTop']), 8):>10} "
                f"{fmt(mp.mpf(row['relativeBErrorMax']), 8):>11} "
                f"{fmt(mp.mpf(row['topNormalizedTErrorMax']), 8):>11}",
                flush=True,
            )
        else:
            print(f"  {grid:4d} {row['sourceSideDim']:3d} active_dim={row['activeDim']}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "point": f(mp.mpf(args.point)),
        "orders": args.orders,
        "basis": args.basis,
        "activeTol": f(mp.mpf(args.active_tol)),
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceGrids": grids,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Weighted source rows approximate the continuum source Hilbert "
            "space.  Stability of active_dim=2, eta, and the perturbative "
            "allowances supports the source-side rank theorem in L2 source "
            "variables."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
