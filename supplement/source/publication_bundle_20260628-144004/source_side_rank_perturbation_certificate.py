#!/usr/bin/env python3
r"""Perturbative certificate for source-side noncollapse.

The analytic rank theorem is

    rank(L E^* | U_delta) = 2,

where ``U_delta`` is the two-dimensional active spectral subspace of
``T = E E^*``.  This script computes the finite source-side rank margin and
the operator-error tolerances which would make the continuum theorem follow by
a Riesz-projector perturbation argument.

Let ``T_N=E_N E_N^*``, ``B_N=L E_N^*``, and let ``Q_N`` be the active projector
onto the top two source eigenvalues.  Put

    m_N = sigma_min(B_N Q_N),
    b_N = ||B_N||,
    g_N = gap from the active pair to the remaining source spectrum.

If the continuum operators satisfy

    ||T-T_N|| < g_N/4,
    ||B-B_N|| + b_N * 4||T-T_N||/g_N < m_N,

then ``B`` has rank two on the continuum active subspace.  The split sufficient
conditions reported here are

    ||B-B_N|| <= m_N/4,
    ||T-T_N||/||T_N|| <= (m_N/b_N) (g_N/||T_N||) / 16.

These are not new numerical evidence for positivity; they are the exact
``finish line'' constants for the source-side noncollapse proof.
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


def op_norm_from_gram(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat * mat.T), eigvals_only=True)
    return mp.sqrt(max([mp.mpf("0"), *vals]))


def singular_values(mat):
    vals = eigvals_gram(mat)
    return [mp.sqrt(max(mp.mpf("0"), val)) for val in vals]


def source_case(args, basis, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol))
    if len(active_idx) != 2:
        raise ValueError(f"basis {basis} has active dimension {len(active_idx)}, expected 2")

    nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs_at_s0)
    source_on_scaled = source_rows * base["scaledModes"]
    source_side = sym(source_on_scaled * source_on_scaled.T)
    source_vals, source_vecs = mp.eigsy(source_side, eigvals_only=False)
    active_left_idx = [i for i, val in enumerate(source_vals) if val > active_floor]
    if len(active_left_idx) != 2:
        raise ValueError(
            f"basis {basis} source-side active dimension {len(active_left_idx)}, expected 2"
        )
    active_left = columns(source_vecs, active_left_idx)

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(point_args, base, base["scaledModes"], mp.mpf(args.point), max(args.orders))
    lmat = submatrix_rows(full_response, args.orders)
    bmat = lmat * source_on_scaled.T
    compressed = bmat * active_left

    raw_svals = singular_values(compressed)
    raw_min = min(raw_svals)
    raw_max = max(raw_svals)
    b_norm = op_norm_from_gram(bmat)
    normalized_response, response_norms = normalized_columns(compressed)
    normalized_svals = singular_values(normalized_response)

    top = max(abs(val) for val in source_vals)
    active_min_idx = active_left_idx[0]
    next_below = source_vals[active_min_idx - 1] if active_min_idx > 0 else mp.mpf("0")
    active_min = source_vals[active_min_idx]
    spectral_gap = active_min - next_below
    gap_normed = spectral_gap / top if top else mp.mpf("0")

    relative_rank_margin = raw_min / b_norm if b_norm else mp.mpf("0")
    split_b_absolute = raw_min / 4
    split_b_relative = relative_rank_margin / 4
    split_t_absolute = raw_min * spectral_gap / (16 * b_norm) if b_norm else mp.mpf("0")
    split_t_top_relative = relative_rank_margin * gap_normed / 16

    return {
        "basis": basis,
        "sourceSideDim": source_side.rows,
        "sourceNodes": [f(x) for x in nodes],
        "activeFloor": f(active_floor),
        "sourceTop": f(top),
        "activeSourceEigenvalues": [f(source_vals[i]) for i in active_left_idx],
        "nextBelowActiveEigenvalue": f(next_below),
        "spectralGapToComplement": f(spectral_gap),
        "spectralGapToComplementFracOfTop": f(gap_normed),
        "rawCompressedSingularValues": [f(x) for x in raw_svals],
        "rawRankMargin": f(raw_min),
        "rawCompressedOpNorm": f(raw_max),
        "responseOperatorNorm": f(b_norm),
        "relativeRankMargin": f(relative_rank_margin),
        "responseColumnNorms": [f(x) for x in response_norms],
        "normalizedResponseSingularValues": [f(x) for x in normalized_svals],
        "normalizedRankMargin": f(min(normalized_svals)),
        "splitSufficientConditions": {
            "absoluteBErrorMax": f(split_b_absolute),
            "relativeBErrorMax": f(split_b_relative),
            "absoluteTErrorMax": f(split_t_absolute),
            "topNormalizedTErrorMax": f(split_t_top_relative),
            "requiresTErrorBelowGapQuarter": f(spectral_gap / 4),
            "requiresTErrorBelowGapQuarterFracOfTop": f(gap_normed / 4),
        },
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
    parser.add_argument("--json-out", default="source_side_rank_perturbation_certificate.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    print(f"Source-side rank perturbation point={args.point} orders={args.orders}", flush=True)
    print("  basis norm_eta rel_margin gap/top allow_B_rel allow_T/top", flush=True)
    rows = []
    for basis in bases:
        row = source_case(args, basis, e_derivs_at_s0)
        rows.append(row)
        cond = row["splitSufficientConditions"]
        print(
            f"  {basis:5d} "
            f"{fmt(mp.mpf(row['normalizedRankMargin']), 8):>10} "
            f"{fmt(mp.mpf(row['relativeRankMargin']), 8):>10} "
            f"{fmt(mp.mpf(row['spectralGapToComplementFracOfTop']), 8):>10} "
            f"{fmt(mp.mpf(cond['relativeBErrorMax']), 8):>11} "
            f"{fmt(mp.mpf(cond['topNormalizedTErrorMax']), 8):>11}",
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
        "rows": rows,
        "theorem": (
            "If ||T-T_N||<gap_N/4 and ||B-B_N|| + ||B_N||*4||T-T_N||/gap_N "
            "< sigma_min(B_N Q_N), then rank(B|U_delta)=2.  The split "
            "sufficient conditions in each row reserve one quarter of the rank "
            "margin for B-error and one quarter for active-projector error."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
