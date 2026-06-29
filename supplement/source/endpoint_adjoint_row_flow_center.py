#!/usr/bin/env python3
r"""Endpoint adjoint row-flow center for the full-rank shortcut.

The endpoint map can be written without propagating the full solution basis:

    M = R_R(s0) - R_L(s0),

where endpoint covectors are transported to ``s0`` by the adjoint row flow

    r'(s) = -r(s)A(s).

This script computes those transported endpoint rows directly.  It then tests
full active row rank by scanning the stable 2x2 minors of M, rather than by
bounding the entire unstable 8-column map in Frobenius norm.

This is a numerical center/refinement test, not an interval proof.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    boundary_row_on_polys,
    companion_matrix,
    mat_frobenius,
    matrix_rank_gram,
    row_norm,
    status,
)
from endpoint_flow_chebyshev_center import (  # noqa: E402
    cheb_lobatto_diff_matrix,
    cheb_lobatto_nodes,
    exact_derivative_at,
)
from global_trace_observability_gap import f, fmt  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def node_key(value: mp.mpf) -> str:
    return mp.nstr(value, 80)


def vec_dot(a, b) -> mp.mpf:
    return mp.fsum(a[i] * b[i] for i in range(len(a)))


def flip_derivs(e_derivs):
    return [[-value for value in row] for row in e_derivs]


def aligned_derivative_cache(args, nodes: list[mp.mpf], max_q: int):
    unique = sorted({node_key(node): node for node in nodes}.values())
    cache = {}
    previous = None
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    rows = []
    for idx, node in enumerate(unique):
        vals, e_derivs, lam_derivs = exact_derivative_at(args, node, max_q)
        if previous is not None:
            cos = vec_dot(previous, e_derivs[0])
            if cos < 0:
                e_derivs = flip_derivs(e_derivs)
                cos = -cos
            min_cos = min(min_cos, abs(cos))
        previous = list(e_derivs[0])
        gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
        min_gap = min(min_gap, gap)
        cache[node_key(node)] = e_derivs
        rows.append(
            {
                "index": idx,
                "s": f(node),
                "lambda0": f(vals[0]),
                "lambda1": f(vals[1]),
                "gap": f(gap),
                "e8": f(e_derivs[0][-1]),
            }
        )
    return cache, rows, min_gap, min_cos


def nearest_derivs(cache: dict[str, list[list[mp.mpf]]], value: mp.mpf):
    key = min(cache, key=lambda candidate: abs(mp.mpf(candidate) - value))
    return cache[key]


def scaling_diag(nodes: list[mp.mpf], mode: str, dim: int) -> mp.matrix:
    scale = mp.matrix(dim, dim)
    if mode == "raw":
        for j in range(dim):
            scale[j, j] = 1
        return scale
    if mode != "taylor":
        raise ValueError(mode)
    radius = (max(nodes) - min(nodes)) / 2
    radius = max(abs(radius), mp.mpf("1e-30"))
    for j in range(dim):
        scale[j, j] = radius**j / mp.factorial(j)
    return scale


def active_endpoint_covectors(polys, active_modes, point, e_derivs, jet_order: int):
    dim = jet_order - 1
    out = mp.matrix(active_modes.cols, dim)
    for col in range(dim):
        unit = [mp.mpf("0") for _ in range(dim)]
        unit[col] = mp.mpf("1")
        brow = trace_green_concomitant_row(e_derivs, unit, jet_order)
        row = boundary_row_on_polys(polys, point, brow)
        active = row * active_modes
        for active_idx in range(active_modes.cols):
            out[active_idx, col] = active[0, active_idx]
    return out


def row_flow_collocation(
    nodes: list[mp.mpf],
    derivs,
    boundary_index: int,
    target_index: int,
    boundary_rows: mp.matrix,
    scaling: str,
) -> mp.matrix:
    order = len(nodes)
    dim = boundary_rows.cols
    row_count = boundary_rows.rows
    left = min(nodes)
    right = max(nodes)
    dmat = cheb_lobatto_diff_matrix(order, left, right)
    scale = scaling_diag(nodes, scaling, dim)
    scale_inv = scale ** -1
    amats = [
        scale * companion_matrix(derivs[j], dim) * scale_inv
        for j in range(order)
    ]
    boundary_scaled = boundary_rows * scale_inv
    size = order * dim
    system = mp.matrix(size, size)

    def idx(node: int, comp: int) -> int:
        return node * dim + comp

    for j in range(order):
        for m in range(dim):
            row = idx(j, m)
            if j == boundary_index:
                system[row, idx(j, m)] = 1
                continue
            for k in range(order):
                system[row, idx(k, m)] += dmat[j, k]
            amat = amats[j]
            for ell in range(dim):
                system[row, idx(j, ell)] += amat[ell, m]

    transported_scaled = mp.matrix(row_count, dim)
    for r in range(row_count):
        rhs = mp.matrix(size, 1)
        for m in range(dim):
            rhs[idx(boundary_index, m)] = boundary_scaled[r, m]
        sol = mp.lu_solve(system, rhs)
        for m in range(dim):
            transported_scaled[r, m] = sol[idx(target_index, m)]
    return transported_scaled * scale


def matrix_to_lists(mat: mp.matrix) -> list[list[float]]:
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def scan_minors(mat: mp.matrix) -> dict:
    minors = []
    best = None
    for i in range(mat.cols):
        for j in range(i + 1, mat.cols):
            det = mat[0, i] * mat[1, j] - mat[0, j] * mat[1, i]
            item = {
                "columns": [i, j],
                "determinant": f(det),
                "absDeterminant": f(abs(det)),
            }
            minors.append(item)
            if best is None or abs(det) > mp.mpf(str(best["absDeterminant"])):
                best = item
    return {
        "bestMinor": best,
        "topMinors": sorted(minors, key=lambda row: row["absDeterminant"], reverse=True)[:8],
    }


def build_center(args, order: int) -> dict:
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    max_q = max(args.max_trace_q, args.jet_order - 1)
    left_nodes = cheb_lobatto_nodes(order, left, s0)
    right_nodes = cheb_lobatto_nodes(order, s0, right)
    all_nodes = left_nodes + right_nodes
    cache, derivative_rows, min_gap, min_cos = aligned_derivative_cache(
        args,
        all_nodes,
        max_q,
    )
    left_derivs = [cache[node_key(node)] for node in left_nodes]
    right_derivs = [cache[node_key(node)] for node in right_nodes]
    e_s0 = nearest_derivs(cache, s0)
    e_left = nearest_derivs(cache, left)
    e_right = nearest_derivs(cache, right)

    base, active_modes, active_idx = active_setup(args, e_s0)
    c_left = active_endpoint_covectors(
        base["polys"],
        active_modes,
        left,
        e_left,
        args.jet_order,
    )
    c_right = active_endpoint_covectors(
        base["polys"],
        active_modes,
        right,
        e_right,
        args.jet_order,
    )
    # Lobatto nodes are high-to-low.  On [left,s0], left endpoint is order-1,
    # s0 is 0.  On [s0,right], right endpoint is 0, s0 is order-1.
    transported_left = row_flow_collocation(
        left_nodes,
        left_derivs,
        order - 1,
        0,
        c_left,
        args.state_scaling,
    )
    transported_right = row_flow_collocation(
        right_nodes,
        right_derivs,
        0,
        order - 1,
        c_right,
        args.state_scaling,
    )
    endpoint_mat = transported_right - transported_left
    rank, eigs, rank_threshold = matrix_rank_gram(
        endpoint_mat * endpoint_mat.T,
        mp.mpf(args.trace_tol),
    )
    sigmas = [mp.sqrt(max(mp.mpf("0"), value)) for value in eigs]
    minor_data = scan_minors(endpoint_mat)
    return {
        "chebyshevOrder": order,
        "stateScaling": args.state_scaling,
        "activeDim": len(active_idx),
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "endpointMapRank": rank,
        "endpointMapGramEigenvalues": [f(value) for value in eigs],
        "endpointMapRankThreshold": f(rank_threshold),
        "sigmaMin": f(min(sigmas)),
        "sigmaMax": f(max(sigmas)),
        "endpointMap": matrix_to_lists(endpoint_mat),
        "endpointMapFrobenius": f(mat_frobenius(endpoint_mat)),
        "transportedLeftRowsFrobenius": f(mat_frobenius(transported_left)),
        "transportedRightRowsFrobenius": f(mat_frobenius(transported_right)),
        "bestMinor": minor_data["bestMinor"],
        "topMinors": minor_data["topMinors"],
        "derivativeRows": derivative_rows,
        "statuses": {
            "activeEndpointMapFullRowRank": status(
                "active endpoint map full row rank",
                rank == len(active_idx),
                "Adjoint row-flow endpoint center has full sampled active row rank.",
            ),
        },
    }


def mp_matrix(rows: list[list[float]]) -> mp.matrix:
    return mp.matrix([[mp.mpf(str(value)) for value in row] for row in rows])


def compare_centers(left: dict, right: dict, rank_margin: mp.mpf) -> dict:
    a = mp_matrix(left["endpointMap"])
    b = mp_matrix(right["endpointMap"])
    diff = b - a
    frob = mat_frobenius(diff)
    left_minor = left["bestMinor"]
    right_minor = right["bestMinor"]
    same_minor = left_minor["columns"] == right_minor["columns"]
    right_det = mp.mpf(str(right_minor["determinant"]))
    left_det_same = None
    det_diff = None
    rel_det_diff = None
    if same_minor:
        left_det_same = mp.mpf(str(left_minor["determinant"]))
        det_diff = abs(right_det - left_det_same)
        rel_det_diff = det_diff / max(mp.mpf("1"), abs(right_det))
    return {
        "leftOrder": left["chebyshevOrder"],
        "rightOrder": right["chebyshevOrder"],
        "frobeniusDifference": f(frob),
        "rankMargin": f(rank_margin),
        "differenceOverMargin": f(frob / rank_margin),
        "passesRankMargin": bool(frob < rank_margin),
        "leftBestMinor": left_minor,
        "rightBestMinor": right_minor,
        "sameBestMinor": same_minor,
        "bestMinorDetAbsDiff": f(det_diff) if det_diff is not None else None,
        "bestMinorRelativeDiff": f(rel_det_diff) if rel_det_diff is not None else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", default="5 7 9")
    parser.add_argument("--rank-margin", default="2.77565955249e6")
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--cache-dir", default=".endpoint_flow_cache")
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
    parser.add_argument("--json-out", default="endpoint_adjoint_row_flow_center.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    orders = parse_ints(args.orders)
    rank_margin = mp.mpf(args.rank_margin)
    centers = []
    for order in orders:
        print(f"Building endpoint adjoint row-flow center order={order}", flush=True)
        center = build_center(args, order)
        centers.append(center)
        print(
            f"  order={order} rank={center['endpointMapRank']} "
            f"||M||_F={fmt(mp.mpf(center['endpointMapFrobenius']), 12)} "
            f"sigma_min={fmt(mp.mpf(center['sigmaMin']), 12)} "
            f"best_minor={center['bestMinor']['absDeterminant']:.12e}",
            flush=True,
        )
    comparisons = [
        compare_centers(left, right, rank_margin)
        for left, right in zip(centers[:-1], centers[1:])
    ]
    frob_stable = bool(comparisons and all(row["passesRankMargin"] for row in comparisons))
    minor_stable = bool(
        comparisons
        and all(
            row["sameBestMinor"]
            and row["bestMinorRelativeDiff"] is not None
            and float(row["bestMinorRelativeDiff"]) < 0.25
            for row in comparisons
        )
    )
    data = {
        "theoremName": "endpoint adjoint row-flow center",
        "interval": [f(mp.mpf(args.constraint_min)), f(mp.mpf(args.constraint_max))],
        "s0": f(mp.mpf(args.s0)),
        "orders": orders,
        "stateScaling": args.state_scaling,
        "rankMargin": f(rank_margin),
        "centers": centers,
        "comparisons": comparisons,
        "frobeniusStabilityStatus": status(
            "row-flow Frobenius stability",
            frob_stable,
            "Closed if successive row-flow endpoint maps differ below the old rank margin.",
        ),
        "minorStabilityStatus": status(
            "row-flow best-minor stability",
            minor_stable,
            (
                "Closed if the same best 2x2 minor persists and its relative "
                "change is below 25 percent across refinements."
            ),
        ),
        "nextTarget": (
            "If a stable minor emerges, build a ball certificate for that "
            "minor instead of a Frobenius certificate for all of M."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint adjoint row-flow center")
    for row in comparisons:
        print(
            f"  {row['leftOrder']}->{row['rightOrder']} "
            f"diff={row['frobeniusDifference']:.12e} "
            f"diff/margin={row['differenceOverMargin']:.12e} "
            f"same_minor={row['sameBestMinor']} "
            f"minor_rel={row['bestMinorRelativeDiff']}",
            flush=True,
        )
    print(f"  Frobenius stable: {frob_stable}")
    print(f"  minor stable: {minor_stable}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
