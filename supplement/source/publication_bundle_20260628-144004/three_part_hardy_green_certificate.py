#!/usr/bin/env python3
r"""Three-part finite certificate for the Hardy/Green proof program.

The continuum target has three pieces:

1. finite low/mid spectral Schur comparison;
2. high-frequency Hardy tail estimate;
3. endpoint trace control preventing low-mode singularity.

This script measures the corresponding finite quantities for

    A = K|ker(R),        B = K(ker(R), U).

For an A-eigenbasis e_j, it reports:

    C_head(M) = || A_{<M}^{-1/2} P_{<M} B ||^2,
    C_tail(M) = || A_{>=M}^{-1/2} P_{>=M} B ||^2,

and the single-mode source ratios

    || <B ., e_j> ||^2 / lambda_j.

It also evaluates the low eigenmodes against a denser Lambda_a grid.  The
sampled constraints force R e_j = 0 at the sample centers; dense leakage shows
how close that finite condition is to the intended continuum moving trace.
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


def make_case(base, basis, constraints=None):
    if constraints is None:
        if base.constraint_rule == "half":
            constraints = max(2, basis // 2)
        elif base.constraint_rule == "basis":
            constraints = basis
        else:
            constraints = base.constraints
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
    return columns(rvecs, n_idx), columns(rvecs, u_idx), rvals


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


def compute_case(base, basis):
    args = make_case(base, basis)
    length = mp.mpf(args.L)
    K, polys = gram_matrix(args, mp.mpf(args.omega), length)
    _, R = trace_matrix(polys, args)
    N, U, _ = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    B = N.T * K * U

    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    full_modes = N * avecs
    beig = avecs.T * B if B.cols else mp.matrix(len(keep), 0)

    gamma2 = window_constant(beig, avals, 0, len(avals))
    split_rows = []
    cutoffs = sorted(set([m for m in parse_ints(base.cutoffs) if 0 <= m <= len(avals)] + list(range(len(avals) + 1))))
    for m in cutoffs:
        head = window_constant(beig, avals, 0, m)
        tail = window_constant(beig, avals, m, len(avals))
        split_rows.append(
            {
                "cutoff": m,
                "head": head,
                "tail": tail,
                "headFrac": head / gamma2 if gamma2 else mp.mpf("0"),
                "tailFrac": tail / gamma2 if gamma2 else mp.mpf("0"),
                "additiveOverhead": (head + tail) / gamma2 if gamma2 else mp.mpf("0"),
            }
        )

    low_count = min(base.low_modes, len(avals))
    dense_args = make_case(base, basis, constraints=base.dense_constraints)
    _, Rdense = trace_matrix(polys, dense_args)
    dense_trace = Rdense * full_modes
    sampled_trace = R * full_modes
    low_rows = []
    for j in range(low_count):
        source2 = row_norm2(beig, j)
        single = source2 / avals[j]
        dense_max, dense_rms = matrix_abs_stats(dense_trace, j)
        sampled_max, sampled_rms = matrix_abs_stats(sampled_trace, j)
        low_rows.append(
            {
                "mode": j,
                "lambda": avals[j],
                "sourceNorm2": source2,
                "singleConstant": single,
                "singleFrac": single / gamma2 if gamma2 else mp.mpf("0"),
                "sourceOverSqrtLambda": mp.sqrt(source2 / avals[j]) if avals[j] > 0 else mp.mpf("0"),
                "sampledTraceMax": sampled_max,
                "sampledTraceRms": sampled_rms,
                "denseTraceMax": dense_max,
                "denseTraceRms": dense_rms,
            }
        )

    tail_thresholds = []
    for threshold in [mp.mpf("0.5"), mp.mpf("0.25"), mp.mpf("0.1"), mp.mpf("0.05")]:
        hit = None
        for row in split_rows:
            if row["tailFrac"] <= threshold:
                hit = row["cutoff"]
                break
        tail_thresholds.append({"threshold": threshold, "cutoff": hit})

    return {
        "basis": basis,
        "constraints": args.constraints,
        "denseConstraints": base.dense_constraints,
        "rank": U.cols,
        "nullity": N.cols,
        "positiveModes": len(avals),
        "gamma2": gamma2,
        "lambdaMin": avals[0] if avals else mp.mpf("0"),
        "lambdaMax": avals[-1] if avals else mp.mpf("0"),
        "splits": split_rows,
        "lowModes": low_rows,
        "tailThresholds": tail_thresholds,
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
        "splits": [
            {
                "cutoff": item["cutoff"],
                "head": f(item["head"]),
                "tail": f(item["tail"]),
                "headFrac": f(item["headFrac"]),
                "tailFrac": f(item["tailFrac"]),
                "additiveOverhead": f(item["additiveOverhead"]),
            }
            for item in row["splits"]
        ],
        "lowModes": [
            {
                "mode": item["mode"],
                "lambda": f(item["lambda"]),
                "sourceNorm2": f(item["sourceNorm2"]),
                "singleConstant": f(item["singleConstant"]),
                "singleFrac": f(item["singleFrac"]),
                "sourceOverSqrtLambda": f(item["sourceOverSqrtLambda"]),
                "sampledTraceMax": f(item["sampledTraceMax"]),
                "sampledTraceRms": f(item["sampledTraceRms"]),
                "denseTraceMax": f(item["denseTraceMax"]),
                "denseTraceRms": f(item["denseTraceRms"]),
            }
            for item in row["lowModes"]
        ],
        "tailThresholds": [
            {
                "threshold": f(item["threshold"]),
                "cutoff": item["cutoff"],
            }
            for item in row["tailThresholds"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="10 12 14 16")
    parser.add_argument("--constraint-rule", choices=["half", "fixed", "basis"], default="half")
    parser.add_argument("--constraints", type=int, default=8)
    parser.add_argument("--dense-constraints", type=int, default=24)
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
    parser.add_argument("--cutoffs", default="4 6 7")
    parser.add_argument("--low-modes", type=int, default=6)
    parser.add_argument("--json-out", default="three_part_hardy_green_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print(
        f"Three-part Hardy/Green certificate model={args.model} kind={args.kind} "
        f"omega={args.omega} L={args.L}",
        flush=True,
    )
    print(
        "  basis cons gamma2      tail<=.5 tail<=.25 tail<=.1 tail<=.05 "
        "low0_frac dense0",
        flush=True,
    )
    for basis in parse_ints(args.basis):
        row = compute_case(args, basis)
        rows.append(serial(row))
        thresholds = {item["threshold"]: item["cutoff"] for item in row["tailThresholds"]}
        low0 = row["lowModes"][0] if row["lowModes"] else None
        print(
            f"  {basis:5d} {row['constraints']:4d} {fmt(row['gamma2']):>10} "
            f"{str(thresholds[mp.mpf('0.5')]):>8} "
            f"{str(thresholds[mp.mpf('0.25')]):>9} "
            f"{str(thresholds[mp.mpf('0.1')]):>8} "
            f"{str(thresholds[mp.mpf('0.05')]):>9} "
            f"{fmt(low0['singleFrac'] if low0 else 0):>9} "
            f"{fmt(low0['denseTraceMax'] if low0 else 0):>10}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
