#!/usr/bin/env python3
r"""Off-base derivative-rank certificate for active trace unique continuation.

If a trace response F_v(a)=Lambda_a(v) vanishes on an open interval, then all
off-base derivatives D_a^q F_v(a0) vanish.  Therefore, on a finite active
space, the implication

    F_v = 0 on an interval  =>  v = 0

is certified by full column rank of the derivative matrix

    M_{q,k}=D_a^q Lambda_a(v_k)|_{a=a0}

at one off-base point.  This is the local analytic/unique-continuation version
of the exterior determinant certificate.
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

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import (  # noqa: E402
    exact_trace_derivatives,
    tower_row_on_polys,
)
from quotient_factorization_mp import columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def parse_mpf_list(text: str) -> list[mp.mpf]:
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def normalized_columns(mat):
    out = mp.matrix(mat.rows, mat.cols)
    norms = []
    for j in range(mat.cols):
        norm = mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows)))
        norms.append(norm)
        scale = 1 / norm if norm else mp.mpf("0")
        for i in range(mat.rows):
            out[i, j] = scale * mat[i, j]
    return out, norms


def submatrix_rows(mat, rows):
    out = mp.matrix(len(rows), mat.cols)
    for i, row in enumerate(rows):
        for j in range(mat.cols):
            out[i, j] = mat[row, j]
    return out


def best_row_minor(mat):
    dim = mat.cols
    best = None
    for rows in itertools.combinations(range(mat.rows), dim):
        det = mp.det(submatrix_rows(mat, rows))
        item = {
            "derivativeOrders": list(rows),
            "determinant": f(det),
            "absDeterminant": f(abs(det)),
        }
        if best is None or abs(det) > mp.mpf(best["absDeterminant"]):
            best = item
    return best


def eigvals_gram(mat):
    if mat.cols == 0:
        return []
    vals = mp.eigsy((mat.T * mat + mat.T * mat) / 2, eigvals_only=True)
    return vals


def derivative_matrix(args, base, active_modes, point: mp.mpf, max_deriv: int):
    point_args = SimpleNamespace(**vars(args))
    point_args.s0 = str(point)
    _vals, e_derivs, _lam_derivs = exact_trace_derivatives(point_args, max_deriv)
    mat = mp.matrix(max_deriv + 1, active_modes.cols)
    for q in range(max_deriv + 1):
        row = tower_row_on_polys(base["polys"], point, e_derivs, q, args.jet_order)
        active_row = row * active_modes
        for j in range(active_modes.cols):
            mat[q, j] = active_row[0, j]
    return mat


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--points", default="0.085625 0.348125 0.545")
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
    parser.add_argument("--json-out", default="trace_active_derivative_rank.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    points = parse_mpf_list(args.points)
    max_q = max(args.max_trace_q, args.basis - args.local_constraint_offset - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs)
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol))
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis

    print(
        f"Trace active derivative-rank basis={args.basis} "
        f"active_dim={len(active_idx)} max_deriv={args.max_deriv}",
        flush=True,
    )
    print("  point rank min_sval_norm best_minor orders", flush=True)

    rows = []
    for point in points:
        mat = derivative_matrix(args, base, active_modes, point, args.max_deriv)
        normalized, col_norms = normalized_columns(mat)
        vals_raw = eigvals_gram(mat)
        vals_norm = eigvals_gram(normalized)
        rank = sum(1 for val in vals_norm if val > mp.mpf(args.rank_tol))
        best = best_row_minor(normalized)
        min_sval_norm = mp.sqrt(max(mp.mpf("0"), vals_norm[0])) if vals_norm else mp.mpf("0")
        rows.append(
            {
                "point": f(point),
                "activeDim": len(active_idx),
                "maxDeriv": args.max_deriv,
                "rankNormalized": rank,
                "columnNorms": [f(x) for x in col_norms],
                "rawGramEigenvalues": [f(x) for x in vals_raw],
                "normalizedGramEigenvalues": [f(x) for x in vals_norm],
                "minSingularValueNormalized": f(min_sval_norm),
                "bestNormalizedDerivativeMinor": best,
            }
        )
        print(
            f"  {fmt(point, 8):>8} {rank:4d} {fmt(min_sval_norm, 8):>13} "
            f"{fmt(mp.mpf(best['absDeterminant']), 8):>11} {best['derivativeOrders']}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Full normalized derivative rank at an off-base point proves the "
            "finite active unique-continuation implication: if Lambda_a(v) "
            "vanishes on an interval, then all derivatives at that point "
            "vanish, hence v=0 on the active block."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
