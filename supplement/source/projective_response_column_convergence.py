#!/usr/bin/env python3
r"""Direct projective response-column convergence diagnostic.

This is the finite model for proving

    \hat M_N -> \hat M

without using the large high-block row norm bound.  For a fixed off-base tuple
``a0`` and derivative orders, it computes the normalized response columns

    \hat x_{j,N}=L v_{j,N}/||L v_{j,N}||

in the common response space R^d, then compares them across Galerkin bases up
to the unavoidable sign and permutation choices of simple eigenvectors.
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

from determinant_gap_bound_diagnostic import (  # noqa: E402
    eigvals_gram,
    parse_orders,
    submatrix_rows,
)
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from trace_active_derivative_rank import derivative_matrix, normalized_columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def matrix_to_lists(mat):
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


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
        signed_dots = []
        for j, k in enumerate(perm):
            dist, sign, c = signed_distance(col(left, j), col(right, k))
            distances.append(dist)
            signs.append(sign)
            signed_dots.append(c)
        max_dist = max(distances) if distances else mp.mpf("0")
        frob = mp.sqrt(mp.fsum(d * d for d in distances))
        item = {
            "permutation": list(perm),
            "signs": signs,
            "signedDots": [f(x) for x in signed_dots],
            "columnDistances": [f(x) for x in distances],
            "maxColumnDistance": f(max_dist),
            "frobeniusDistance": f(frob),
        }
        if best is None or frob < mp.mpf(best["frobeniusDistance"]):
            best = item
    return best


def basis_response(args, basis, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(
        base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol)
    )
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(
        point_args, max(args.orders)
    )
    full_matrix = derivative_matrix(
        point_args, base, active_modes, mp.mpf(args.point), max(args.orders)
    )
    raw = submatrix_rows(full_matrix, args.orders)
    normalized, column_norms = normalized_columns(raw)
    norm_vals = eigvals_gram(normalized)
    eta = mp.sqrt(max(mp.mpf("0"), norm_vals[0])) if norm_vals else mp.mpf("0")
    return {
        "basis": basis,
        "activeDim": len(active_idx),
        "activeFloor": f(active_floor),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "columnNorms": [f(x) for x in column_norms],
        "normalizedSingularValuesSquared": [f(x) for x in norm_vals],
        "normalizedMinSingularValue": f(eta),
        "normalizedResponseColumns": matrix_to_lists(normalized),
        "_normalizedMatrix": normalized,
    }


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
    parser.add_argument("--json-out", default="projective_response_column_convergence.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    print(
        f"Projective response-column convergence point={args.point} "
        f"orders={args.orders} bases={bases}",
        flush=True,
    )
    rows = []
    matrices = {}
    for basis in bases:
        row = basis_response(args, basis, e_derivs_at_s0)
        matrices[basis] = row.pop("_normalizedMatrix")
        rows.append(row)
        print(
            f"  basis {basis:2d}, dim {row['activeDim']}: "
            f"eta={fmt(mp.mpf(row['normalizedMinSingularValue']), 8)} "
            f"cols={row['normalizedResponseColumns']}",
            flush=True,
        )

    comparisons = []
    for left, right in zip(bases[:-1], bases[1:]):
        item = compare_columns(matrices[left], matrices[right])
        item["leftBasis"] = left
        item["rightBasis"] = right
        comparisons.append(item)
        print(
            f"  compare {left}->{right}: max_col={fmt(mp.mpf(item['maxColumnDistance']), 8)} "
            f"frob={fmt(mp.mpf(item['frobeniusDistance']), 8)} "
            f"perm={item['permutation']} signs={item['signs']}",
            flush=True,
        )

    reference = bases[-1]
    reference_comparisons = []
    for basis in bases[:-1]:
        item = compare_columns(matrices[basis], matrices[reference])
        item["leftBasis"] = basis
        item["rightBasis"] = reference
        reference_comparisons.append(item)

    max_consecutive = max(
        [mp.mpf(item["maxColumnDistance"]) for item in comparisons] + [mp.mpf("0")]
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
        "rows": rows,
        "consecutiveComparisons": comparisons,
        "referenceBasis": reference,
        "referenceComparisons": reference_comparisons,
        "maxConsecutiveColumnDistance": f(max_consecutive),
        "interpretation": (
            "Direct projective convergence should make these normalized "
            "response columns Cauchy after sign/permutation alignment.  Large "
            "distances mean the current Galerkin active eigenlines are not "
            "yet a stable sequence in the chosen response coordinates."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
