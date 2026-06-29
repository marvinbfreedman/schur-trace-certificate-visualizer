#!/usr/bin/env python3
r"""Coupled refinement scan for the high-frequency Hardy tail.

The finite proof program now has two moving parts:

  1. sampled endpoint traces must approximate the closed moving trace;
  2. the high spectral tail of A=K|ker(R) must be small for the source B(U).

This script increases the Galerkin basis and the number of trace samples
together.  It reports dense trace leakage and tail constants

    C_tail(M) = ||A_{>=M}^{-1/2} P_{>=M} B||^2

as fractions of the full Douglas constant.
"""

from __future__ import annotations

import argparse
import json
import math
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


def constraints_for(base, basis):
    if base.constraint_rule == "ratio":
        return max(1, min(basis - 1, int(round(base.constraint_ratio * basis))))
    if base.constraint_rule == "target-nullity":
        return max(1, min(basis - 1, basis - base.target_nullity))
    if base.constraint_rule == "offset":
        return max(1, min(basis - 1, basis - base.constraint_offset))
    raise ValueError(base.constraint_rule)


def make_args(base, basis, constraints):
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=base.quad_base + base.quad_step * basis,
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


def matrix_abs_stats(mat, col, limit):
    cols = min(limit, mat.cols)
    best = mp.mpf("0")
    rms_acc = mp.mpf("0")
    count = 0
    for j in range(cols):
        vals = [abs(mat[i, j]) for i in range(mat.rows)]
        if vals:
            best = max(best, max(vals))
            rms_acc += mp.fsum(v * v for v in vals) / len(vals)
            count += 1
    rms = mp.sqrt(rms_acc / count) if count else mp.mpf("0")
    return best, rms


def first_tail_cutoff(tails, threshold):
    for row in tails:
        if row["fraction"] <= threshold:
            return row["start"]
    return None


def compute_case(base, basis):
    constraints = constraints_for(base, basis)
    args = make_args(base, basis, constraints)
    length = mp.mpf(args.L)
    K, polys = gram_matrix(args, mp.mpf(args.omega), length)
    _, R = trace_matrix(polys, args)
    N, U, rank, nullity = split_spaces(R, args.rank_tol)
    if N.cols == 0 or U.cols == 0:
        return None

    A = N.T * K * N
    B = N.T * K * U
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    if not keep:
        return None

    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    modes = N * avecs
    beig = avecs.T * B
    gamma2 = window_constant(beig, avals, 0, len(avals))
    tails = []
    for start in range(len(avals) + 1):
        value = window_constant(beig, avals, start, len(avals))
        tails.append(
            {
                "start": start,
                "constant": value,
                "fraction": value / gamma2 if gamma2 else mp.mpf("0"),
            }
        )

    dense_constraints = base.dense_constraints or max(constraints + 6, int(math.ceil(base.dense_ratio * basis)))
    dense_args = make_args(base, basis, dense_constraints)
    _, dense_R = trace_matrix(polys, dense_args)
    dense_trace = dense_R * modes
    dense_max, dense_rms = matrix_abs_stats(dense_trace, 0, base.low_modes)

    mode_rows = []
    for j in range(len(avals)):
        source2 = row_norm2(beig, j)
        single = source2 / avals[j]
        mode_rows.append(
            {
                "mode": j,
                "lambda": avals[j],
                "singleConstant": single,
                "singleFrac": single / gamma2 if gamma2 else mp.mpf("0"),
            }
        )
    top_modes = sorted(mode_rows, key=lambda row: row["singleFrac"], reverse=True)[: base.top]

    return {
        "basis": basis,
        "constraints": constraints,
        "denseConstraints": dense_constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "gamma2": gamma2,
        "lambdaMin": avals[0],
        "lambdaMax": avals[-1],
        "denseTraceMax": dense_max,
        "denseTraceRms": dense_rms,
        "tailCutoff10": first_tail_cutoff(tails, mp.mpf("0.1")),
        "tailCutoff05": first_tail_cutoff(tails, mp.mpf("0.05")),
        "tails": tails,
        "topModes": top_modes,
    }


def serial(row):
    return {
        "basis": row["basis"],
        "constraints": row["constraints"],
        "denseConstraints": row["denseConstraints"],
        "rank": row["rank"],
        "nullity": row["nullity"],
        "positiveModes": row["positiveModes"],
        "gamma2": f(row["gamma2"]),
        "lambdaMin": f(row["lambdaMin"]),
        "lambdaMax": f(row["lambdaMax"]),
        "denseTraceMax": f(row["denseTraceMax"]),
        "denseTraceRms": f(row["denseTraceRms"]),
        "tailCutoff10": row["tailCutoff10"],
        "tailCutoff05": row["tailCutoff05"],
        "tails": [
            {
                "start": item["start"],
                "constant": f(item["constant"]),
                "fraction": f(item["fraction"]),
            }
            for item in row["tails"]
        ],
        "topModes": [
            {
                "mode": item["mode"],
                "lambda": f(item["lambda"]),
                "singleConstant": f(item["singleConstant"]),
                "singleFrac": f(item["singleFrac"]),
            }
            for item in row["topModes"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="12 14 16 18")
    parser.add_argument("--constraint-rule", choices=["ratio", "target-nullity", "offset"], default="ratio")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--target-nullity", type=int, default=6)
    parser.add_argument("--constraint-offset", type=int, default=6)
    parser.add_argument("--dense-constraints", type=int, default=0)
    parser.add_argument("--dense-ratio", type=float, default=1.6)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--quad-base", type=int, default=6)
    parser.add_argument("--quad-step", type=int, default=1)
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
    parser.add_argument("--json-out", default="high_frequency_tail_refinement.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print(
        f"High-frequency tail refinement model={args.model} kind={args.kind} "
        f"rule={args.constraint_rule}",
        flush=True,
    )
    print(
        "  basis cons null gamma2      leak       tail10 tail05 top_modes",
        flush=True,
    )
    for basis in parse_ints(args.basis):
        row = compute_case(args, basis)
        if row is None:
            continue
        rows.append(serial(row))
        top = ",".join(str(item["mode"]) for item in row["topModes"])
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['nullity']:4d} "
            f"{fmt(row['gamma2']):>10} {fmt(row['denseTraceMax']):>10} "
            f"{str(row['tailCutoff10']):>6} {str(row['tailCutoff05']):>6} {top}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
