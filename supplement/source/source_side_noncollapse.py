#!/usr/bin/env python3
r"""Source-side noncollapse diagnostic.

For S=E^*E, any eigenvector with positive eigenvalue is in Range(E^*).
Therefore response noncollapse can be checked on the source side:

    L v_j = mu_j^{-1/2} L E^* u_j,

where u_j is the corresponding eigenvector of EE^*.  This is the right finite
model after the whole-kernel domination test fails.
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
from global_trace_observability_gap import f, fmt, parse_ints, stack_source_rows  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from trace_active_derivative_rank import derivative_matrix, normalized_columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def matrix_to_lists(mat):
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def source_side_case(args, basis, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(
        base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol)
    )
    active_basis = columns(base["sourceVecs"], active_idx)

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(
        point_args,
        base,
        base["scaledModes"],
        mp.mpf(args.point),
        max(args.orders),
    )
    Lmat = submatrix_rows(full_response, args.orders)
    response_active_right = Lmat * active_basis
    normalized_right, right_col_norms = normalized_columns(response_active_right)
    right_svals2 = eigvals_gram(normalized_right)

    _nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs_at_s0)
    source_on_scaled = source_rows * base["scaledModes"]
    source_side_gram = sym(source_on_scaled * source_on_scaled.T)
    source_vals, source_vecs = mp.eigsy(source_side_gram, eigvals_only=False)
    active_left_idx = [i for i, val in enumerate(source_vals) if val > active_floor]
    active_left = columns(source_vecs, active_left_idx)
    source_response = Lmat * source_on_scaled.T * active_left
    normalized_source, source_col_norms = normalized_columns(source_response)
    source_svals2 = eigvals_gram(normalized_source)
    return {
        "basis": basis,
        "sourceSideDim": source_side_gram.rows,
        "activeDim": len(active_idx),
        "activeFloor": f(active_floor),
        "rightActiveEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "sourceSideActiveEigenvalues": [f(source_vals[i]) for i in active_left_idx],
        "rightResponseColumnNorms": [f(x) for x in right_col_norms],
        "sourceResponseColumnNorms": [f(x) for x in source_col_norms],
        "rightNormalizedMinSingularValue": f(mp.sqrt(max(mp.mpf("0"), right_svals2[0])) if right_svals2 else mp.mpf("0")),
        "sourceNormalizedMinSingularValue": f(mp.sqrt(max(mp.mpf("0"), source_svals2[0])) if source_svals2 else mp.mpf("0")),
        "rightNormalizedColumns": matrix_to_lists(normalized_right),
        "sourceNormalizedColumns": matrix_to_lists(normalized_source),
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
    parser.add_argument("--json-out", default="source_side_noncollapse.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    print(f"Source-side noncollapse point={args.point} orders={args.orders} bases={bases}", flush=True)
    print("  basis dim right_eta source_eta source_norms", flush=True)
    rows = []
    for basis in bases:
        row = source_side_case(args, basis, e_derivs_at_s0)
        rows.append(row)
        print(
            f"  {basis:5d} {row['activeDim']:3d} "
            f"{fmt(mp.mpf(row['rightNormalizedMinSingularValue']), 8):>10} "
            f"{fmt(mp.mpf(row['sourceNormalizedMinSingularValue']), 8):>10} "
            f"{row['sourceResponseColumnNorms']}",
            flush=True,
        )

    min_source_eta = min(mp.mpf(row["sourceNormalizedMinSingularValue"]) for row in rows)
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
        "minSourceSideNormalizedSingularValue": f(min_source_eta),
        "rows": rows,
        "interpretation": (
            "Positive source-side response rank proves finite noncollapse for "
            "positive source eigenvectors because v_j = mu_j^{-1/2} E^*u_j."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
