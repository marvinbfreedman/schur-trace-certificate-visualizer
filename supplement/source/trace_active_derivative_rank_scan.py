#!/usr/bin/env python3
r"""Cross-basis scan for the off-base derivative-rank certificate.

For each Galerkin basis, this script builds the active source block and tests
whether the analytic trace-response derivatives

    D_a^q Lambda_a(v_k)|_{a=a0}

have full column rank at off-base points.  Full rank is the finite-dimensional
unique-continuation certificate: if Lambda_a(v) vanishes on an interval
containing a0, all these rows vanish, so v=0 on the active block.
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
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from trace_active_derivative_rank import (  # noqa: E402
    best_row_minor,
    derivative_matrix,
    eigvals_gram,
    normalized_columns,
    parse_mpf_list,
)
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def parse_int_list(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def basis_args(args, basis: int):
    out = SimpleNamespace(**vars(args))
    out.basis = basis
    return out


def scan_basis(args, basis: int, points: list[mp.mpf]) -> dict:
    bargs = basis_args(args, basis)
    max_q = max(bargs.max_trace_q, basis - bargs.local_constraint_offset - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(bargs, max_q)
    base = local_case(bargs, basis, e_derivs)
    active_idx, active_floor = active_indices(
        base["sourceVals"], base["sourceTop"], mp.mpf(bargs.active_tol)
    )
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis

    point_rows = []
    for point in points:
        mat = derivative_matrix(bargs, base, active_modes, point, bargs.max_deriv)
        normalized, col_norms = normalized_columns(mat)
        vals_raw = eigvals_gram(mat)
        vals_norm = eigvals_gram(normalized)
        rank = sum(1 for val in vals_norm if val > mp.mpf(bargs.rank_tol))
        min_sval_norm = mp.sqrt(max(mp.mpf("0"), vals_norm[0])) if vals_norm else mp.mpf("0")
        best = best_row_minor(normalized)
        point_rows.append(
            {
                "point": f(point),
                "activeDim": len(active_idx),
                "rankNormalized": rank,
                "fullRank": rank == len(active_idx),
                "minSingularValueNormalized": f(min_sval_norm),
                "columnNorms": [f(x) for x in col_norms],
                "rawGramEigenvalues": [f(x) for x in vals_raw],
                "normalizedGramEigenvalues": [f(x) for x in vals_norm],
                "bestNormalizedDerivativeMinor": best,
            }
        )

    def point_score(row: dict) -> mp.mpf:
        return mp.mpf(row["minSingularValueNormalized"])

    full_rank_rows = [row for row in point_rows if row["fullRank"]]
    best_point = max(full_rank_rows or point_rows, key=point_score)
    return {
        "basis": basis,
        "activeDim": len(active_idx),
        "activeTol": f(mp.mpf(bargs.active_tol)),
        "activeFloor": f(active_floor),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "points": point_rows,
        "bestPoint": best_point,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--points", default="0.348125 0.545")
    parser.add_argument("--bases", default="18 20 22 24")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
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
    parser.add_argument("--json-out", default="trace_active_derivative_rank_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_int_list(args.bases)
    points = parse_mpf_list(args.points)

    print(
        "Trace active derivative-rank scan "
        f"bases={bases} points={[fmt(p, 8) for p in points]} max_deriv={args.max_deriv}",
        flush=True,
    )
    rows = []
    for basis in bases:
        row = scan_basis(args, basis, points)
        rows.append(row)
        best = row["bestPoint"]
        minor = best["bestNormalizedDerivativeMinor"]
        print(
            "  basis {basis:2d}, dim {dim}: best point {point}, "
            "rank {rank}, min_sval {sval}, minor {minor}, orders {orders}".format(
                basis=basis,
                dim=row["activeDim"],
                point=best["point"],
                rank=best["rankNormalized"],
                sval=fmt(mp.mpf(best["minSingularValueNormalized"]), 8),
                minor=fmt(mp.mpf(minor["absDeterminant"]), 8),
                orders=minor["derivativeOrders"],
            ),
            flush=True,
        )

    min_best_sval = min(mp.mpf(row["bestPoint"]["minSingularValueNormalized"]) for row in rows)
    all_full_rank = all(row["bestPoint"]["fullRank"] for row in rows)
    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "bases": bases,
        "points": [f(p) for p in points],
        "maxDeriv": args.max_deriv,
        "rankTol": f(mp.mpf(args.rank_tol)),
        "minBestSingularValueNormalized": f(min_best_sval),
        "allBestPointsFullRank": all_full_rank,
        "rows": rows,
        "interpretation": (
            "Stable full derivative rank at off-base points is the finite "
            "Galerkin shadow of active unique continuation.  The continuum "
            "proof must show that these confluent derivative matrices converge "
            "on the spectral active spaces in the A-graph norm."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  min best normalized singular value: {fmt(min_best_sval, 8)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
