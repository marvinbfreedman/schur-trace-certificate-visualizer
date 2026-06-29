#!/usr/bin/env python3
r"""Refine sampled endpoint traces at fixed Galerkin basis.

The previous three-part certificate showed that finite sampled constraints can
annihilate the chosen trace rows while still leaking against a denser moving
trace grid.  This script fixes the Galerkin basis and increases the number of
sampled Lambda_a constraints, measuring:

  * the Douglas constant on the quotient split;
  * the high-frequency tail cutoff;
  * source coupling by A-eigenmode;
  * dense moving-trace leakage of the constrained modes.

This is a finite diagnostic for the continuum endpoint trace theorem.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import columns, gram_matrix, max_eig_or_zero, trace_matrix  # noqa: E402


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def f(x):
    return float(x)


def fmt(x, digits=8):
    return mp.nstr(x, digits)


def make_args(base, constraints):
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=base.basis,
        quad=base.quad,
        laguerre=base.laguerre,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
        dps=base.dps,
    )


def split_spaces(R, rank_tol_text):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx), len(u_idx), len(n_idx)


def row_norm2(mat, row_index):
    return mp.fsum(mat[row_index, j] ** 2 for j in range(mat.cols))


def row_outer(row, scale):
    out = mp.matrix(row.cols)
    for i in range(row.cols):
        for j in range(row.cols):
            out[i, j] = scale * row[0, i] * row[0, j]
    return out


def window_constant(beig, avals, start, stop):
    if beig.cols == 0 or start >= stop:
        return mp.mpf("0")
    mat = mp.matrix(beig.cols)
    for j in range(start, stop):
        row = mp.matrix(1, beig.cols)
        for col in range(beig.cols):
            row[0, col] = beig[j, col]
        mat += row_outer(row, 1 / avals[j])
    return max_eig_or_zero(mat)


def matrix_abs_stats(mat, col):
    vals = [abs(mat[i, col]) for i in range(mat.rows)]
    max_abs = max(vals) if vals else mp.mpf("0")
    rms = mp.sqrt(mp.fsum(v * v for v in vals) / len(vals)) if vals else mp.mpf("0")
    return max_abs, rms


def threshold_cutoff(beig, avals, gamma2, threshold):
    if gamma2 == 0:
        return None
    for start in range(len(avals) + 1):
        tail = window_constant(beig, avals, start, len(avals))
        if tail / gamma2 <= threshold:
            return start
    return None


def compute_for_constraints(base, K, polys, dense_R, constraints):
    args = make_args(base, constraints)
    _, R = trace_matrix(polys, args)
    N, U, rank, nullity = split_spaces(R, args.rank_tol)
    if N.cols == 0 or U.cols == 0:
        return {
            "constraints": constraints,
            "rank": rank,
            "nullity": nullity,
            "positiveModes": 0,
            "gamma2": 0.0,
            "tailCutoffs": [],
            "topModes": [],
            "lowModeDenseMax": 0.0,
            "lowModeDenseRms": 0.0,
        }

    A = N.T * K * N
    B = N.T * K * U
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    if not keep:
        return {
            "constraints": constraints,
            "rank": rank,
            "nullity": nullity,
            "positiveModes": 0,
            "gamma2": 0.0,
            "tailCutoffs": [],
            "topModes": [],
            "lowModeDenseMax": 0.0,
            "lowModeDenseRms": 0.0,
        }

    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    full_modes = N * avecs
    beig = avecs.T * B
    dense_trace = dense_R * full_modes
    sampled_trace = R * full_modes
    gamma2 = window_constant(beig, avals, 0, len(avals))

    mode_rows = []
    for j in range(len(avals)):
        source2 = row_norm2(beig, j)
        single = source2 / avals[j]
        dense_max, dense_rms = matrix_abs_stats(dense_trace, j)
        sampled_max, sampled_rms = matrix_abs_stats(sampled_trace, j)
        mode_rows.append(
            {
                "mode": j,
                "lambda": avals[j],
                "singleConstant": single,
                "singleFrac": single / gamma2 if gamma2 else mp.mpf("0"),
                "sourceOverSqrtLambda": mp.sqrt(single),
                "denseTraceMax": dense_max,
                "denseTraceRms": dense_rms,
                "sampledTraceMax": sampled_max,
                "sampledTraceRms": sampled_rms,
            }
        )

    low_count = min(base.low_modes, len(mode_rows))
    low_dense_max = max((mode_rows[j]["denseTraceMax"] for j in range(low_count)), default=mp.mpf("0"))
    low_dense_rms = mp.sqrt(
        mp.fsum(mode_rows[j]["denseTraceRms"] ** 2 for j in range(low_count)) / low_count
    ) if low_count else mp.mpf("0")
    top_modes = sorted(mode_rows, key=lambda row: row["singleFrac"], reverse=True)[: base.top]
    tail_cutoffs = [
        {
            "threshold": threshold,
            "cutoff": threshold_cutoff(beig, avals, gamma2, threshold),
        }
        for threshold in [mp.mpf("0.5"), mp.mpf("0.25"), mp.mpf("0.1"), mp.mpf("0.05")]
    ]

    return {
        "constraints": constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "gamma2": f(gamma2),
        "lambdaMin": f(avals[0]),
        "lambdaMax": f(avals[-1]),
        "tailCutoffs": [
            {"threshold": f(row["threshold"]), "cutoff": row["cutoff"]}
            for row in tail_cutoffs
        ],
        "topModes": [
            {
                "mode": row["mode"],
                "lambda": f(row["lambda"]),
                "singleFrac": f(row["singleFrac"]),
                "sourceOverSqrtLambda": f(row["sourceOverSqrtLambda"]),
                "denseTraceMax": f(row["denseTraceMax"]),
                "denseTraceRms": f(row["denseTraceRms"]),
                "sampledTraceMax": f(row["sampledTraceMax"]),
            }
            for row in top_modes
        ],
        "lowModeDenseMax": f(low_dense_max),
        "lowModeDenseRms": f(low_dense_rms),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=16)
    parser.add_argument("--constraints", default="5 6 7 8 9 10 11 12")
    parser.add_argument("--dense-constraints", type=int, default=32)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--quad", type=int, default=22)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--low-modes", type=int, default=4)
    parser.add_argument("--top", type=int, default=4)
    parser.add_argument("--json-out", default="endpoint_trace_refinement_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    length = mp.mpf(args.L)
    K, polys = gram_matrix(args, mp.mpf(args.omega), length)
    dense_args = make_args(args, args.dense_constraints)
    _, dense_R = trace_matrix(polys, dense_args)

    rows = []
    print(
        f"Endpoint trace refinement scan model={args.model} kind={args.kind} "
        f"basis={args.basis} omega={args.omega} L={args.L}",
        flush=True,
    )
    print(
        "  cons rank null gamma2      tail<=.1 tail<=.05 top modes        "
        "low_dense_max",
        flush=True,
    )
    for constraints in parse_ints(args.constraints):
        row = compute_for_constraints(args, K, polys, dense_R, constraints)
        rows.append(row)
        thresholds = {str(item["threshold"]): item["cutoff"] for item in row["tailCutoffs"]}
        top = ",".join(str(item["mode"]) for item in row["topModes"])
        print(
            f"  {constraints:4d} {row['rank']:4d} {row['nullity']:4d} "
            f"{fmt(row['gamma2']):>10} "
            f"{str(thresholds.get('0.1')):>8} {str(thresholds.get('0.05')):>9} "
            f"{top:>15} {fmt(row['lowModeDenseMax']):>13}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(length),
        "basis": args.basis,
        "denseConstraints": args.dense_constraints,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
