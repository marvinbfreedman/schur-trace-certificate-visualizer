#!/usr/bin/env python3
r"""Stability of the source-side active eigenspace.

After reducing noncollapse to

    rank(L E^* | U_delta)=2,

the natural object is the source-side operator EE^* on the fixed source sample
space.  This script compares the normalized source-side operators, active
source eigenspaces, and normalized response columns across Galerkin bases.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import eigvals_gram, parse_orders, submatrix_rows, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints, stack_source_rows  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from trace_active_derivative_rank import derivative_matrix, normalized_columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def col(mat, j):
    return [mat[i, j] for i in range(mat.rows)]


def dot(u, v):
    return mp.fsum(u[i] * v[i] for i in range(len(u)))


def signed_distance(u, v):
    c = dot(u, v)
    sign = 1 if c >= 0 else -1
    dist = mp.sqrt(mp.fsum((u[i] - sign * v[i]) ** 2 for i in range(len(u))))
    return dist, sign, c


def compare_columns(left, right):
    dim = left.cols
    best = None
    for perm in itertools.permutations(range(dim)):
        distances = []
        signs = []
        dots = []
        for j, k in enumerate(perm):
            dist, sign, c = signed_distance(col(left, j), col(right, k))
            distances.append(dist)
            signs.append(sign)
            dots.append(c)
        frob = mp.sqrt(mp.fsum(d * d for d in distances))
        item = {
            "permutation": list(perm),
            "signs": signs,
            "signedDots": [f(x) for x in dots],
            "columnDistances": [f(x) for x in distances],
            "maxColumnDistance": f(max(distances) if distances else mp.mpf("0")),
            "frobeniusDistance": f(frob),
        }
        if best is None or frob < mp.mpf(best["frobeniusDistance"]):
            best = item
    return best


def eig_op_norm_symmetric(mat):
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def principal_angles(left, right):
    cross = left.T * right
    vals = eigvals_gram(cross)
    svals = [mp.sqrt(max(mp.mpf("0"), v)) for v in vals]
    smin = min(svals) if svals else mp.mpf("0")
    return {
        "singularValues": [f(x) for x in svals],
        "minCosine": f(smin),
        "projectionGap": f(mp.sqrt(max(mp.mpf("0"), 1 - smin**2))),
        "polarAlpha": f(mp.sqrt(max(mp.mpf("0"), 2 - 2 * smin))),
    }


def case(args, basis, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol))

    _nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs_at_s0)
    source_on_scaled = source_rows * base["scaledModes"]
    source_side = sym(source_on_scaled * source_on_scaled.T)
    source_top = max([abs(v) for v in mp.eigsy(source_side, eigvals_only=True)] + [mp.mpf("1")])
    source_normed = source_side / source_top
    vals, vecs = mp.eigsy(source_side, eigvals_only=False)
    active_left_idx = [i for i, val in enumerate(vals) if val > active_floor]
    active_left = columns(vecs, active_left_idx)

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(point_args, base, base["scaledModes"], mp.mpf(args.point), max(args.orders))
    Lmat = submatrix_rows(full_response, args.orders)
    response = Lmat * source_on_scaled.T * active_left
    normalized_response, response_norms = normalized_columns(response)
    response_svals2 = eigvals_gram(normalized_response)
    return {
        "basis": basis,
        "sourceSideDim": source_side.rows,
        "activeDim": len(active_left_idx),
        "activeFloor": f(active_floor),
        "sourceTop": f(source_top),
        "activeSourceEigenvalues": [f(vals[i]) for i in active_left_idx],
        "normalizedActiveEigenvalues": [f(vals[i] / source_top) for i in active_left_idx],
        "responseColumnNorms": [f(x) for x in response_norms],
        "responseNormalizedMinSingularValue": f(mp.sqrt(max(mp.mpf("0"), response_svals2[0])) if response_svals2 else mp.mpf("0")),
        "_sourceNormed": source_normed,
        "_activeLeft": active_left,
        "_normalizedResponse": normalized_response,
    }


def strip_private(row):
    return {key: value for key, value in row.items() if not key.startswith("_")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20 22 24")
    parser.add_argument("--basis", type=int, default=24)
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
    parser.add_argument("--json-out", default="source_side_stability.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    print(f"Source-side stability point={args.point} orders={args.orders} bases={bases}", flush=True)
    rows = []
    by_basis = {}
    for basis in bases:
        row = case(args, basis, e_derivs_at_s0)
        rows.append(row)
        by_basis[basis] = row
        print(
            f"  basis {basis:2d}, dim {row['activeDim']}: "
            f"eta={fmt(mp.mpf(row['responseNormalizedMinSingularValue']), 8)} "
            f"active_normed={row['normalizedActiveEigenvalues']}",
            flush=True,
        )

    comparisons = []
    for left, right in zip(bases[:-1], bases[1:]):
        lrow = by_basis[left]
        rrow = by_basis[right]
        op_diff = eig_op_norm_symmetric(lrow["_sourceNormed"] - rrow["_sourceNormed"])
        angles = principal_angles(lrow["_activeLeft"], rrow["_activeLeft"])
        response = compare_columns(lrow["_normalizedResponse"], rrow["_normalizedResponse"])
        item = {
            "leftBasis": left,
            "rightBasis": right,
            "normalizedOperatorDiff": f(op_diff),
            "activeSubspace": angles,
            "responseColumns": response,
        }
        comparisons.append(item)
        print(
            f"  compare {left}->{right}: op={fmt(op_diff, 8)} "
            f"gap={fmt(mp.mpf(angles['projectionGap']), 8)} "
            f"resp={fmt(mp.mpf(response['maxColumnDistance']), 8)}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "point": f(mp.mpf(args.point)),
        "orders": args.orders,
        "bases": bases,
        "activeTol": f(mp.mpf(args.active_tol)),
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "rows": [strip_private(row) for row in rows],
        "comparisons": comparisons,
        "interpretation": (
            "This is the fixed-source-space stability diagnostic for the "
            "continuum theorem rank(LE^*|U_delta)=2.  The operator difference "
            "and active-subspace gap measure whether EE^* is stabilizing on "
            "the source side."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
