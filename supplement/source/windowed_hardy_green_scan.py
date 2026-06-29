#!/usr/bin/env python3
r"""Scan windowed Hardy/Green constants for the source range B(U).

For the quotient split V = ker(R) \oplus U, write

    A = K|ker(R),        B = K(ker(R), U).

If A e_j = lambda_j e_j, then the spectral-window estimate for a finite
window I is

    sum_{j in I} |<B u, e_j>|^2 / lambda_j <= C_I ||u||_U^2,

where

    C_I = lambda_max(B^* P_I A_I^{-1} P_I B).

This script computes C_I for all contiguous windows.  The full window is the
ordinary Douglas constant; the smaller windows show where the Hardy/Green
proof has to control source projections.
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


def make_case(base, basis):
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
        top=base.top,
    )


def split_spaces(K, R, rank_tol_text):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx)


def row_outer(row, scale):
    out = mp.matrix(row.cols)
    for i in range(row.cols):
        for j in range(row.cols):
            out[i, j] = scale * row[0, i] * row[0, j]
    return out


def window_constant(beig, avals, start, stop):
    if beig.cols == 0:
        return mp.mpf("0")
    mat = mp.matrix(beig.cols)
    for j in range(start, stop):
        row = mp.matrix(1, beig.cols)
        for col in range(beig.cols):
            row[0, col] = beig[j, col]
        mat += row_outer(row, 1 / avals[j])
    return max_eig_or_zero(mat)


def compute_case(args):
    length = mp.mpf(args.L)
    K, polys = gram_matrix(args, mp.mpf(args.omega), length)
    _, R = trace_matrix(polys, args)
    N, U = split_spaces(K, R, args.rank_tol)
    A = N.T * K * N
    B = N.T * K * U
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    beig = avecs.T * B if B.cols else mp.matrix(len(keep), 0)

    windows = []
    cumulative = mp.mpf("0")
    single_constants = []
    tail_constants = []
    for start in range(len(avals)):
        tail_constants.append(window_constant(beig, avals, start, len(avals)))
        for stop in range(start + 1, len(avals) + 1):
            const = window_constant(beig, avals, start, stop)
            windows.append(
                {
                    "start": start,
                    "stop": stop,
                    "width": stop - start,
                    "constant": const,
                    "lambdaStart": avals[start],
                    "lambdaStop": avals[stop - 1],
                }
            )
            if stop == start + 1:
                single_constants.append(const)
            if start == 0 and stop == len(avals):
                cumulative = const

    top_windows = sorted(windows, key=lambda row: row["constant"], reverse=True)[: args.top]
    top_singles = sorted(
        [
            {
                "mode": i,
                "constant": single_constants[i],
                "lambda": avals[i],
            }
            for i in range(len(single_constants))
        ],
        key=lambda row: row["constant"],
        reverse=True,
    )[: args.top]
    tails = [
        {
            "start": i,
            "constant": tail_constants[i],
            "lambdaStart": avals[i],
        }
        for i in range(len(tail_constants))
    ]
    return {
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": U.cols,
        "nullity": N.cols,
        "gamma2Full": cumulative,
        "positiveModes": len(avals),
        "topWindows": top_windows,
        "topSingles": top_singles,
        "tails": tails,
    }


def serial_case(case):
    data = compute_case(case)
    return {
        "basis": data["basis"],
        "constraints": data["constraints"],
        "rank": data["rank"],
        "nullity": data["nullity"],
        "gamma2Full": f(data["gamma2Full"]),
        "positiveModes": data["positiveModes"],
        "topWindows": [
            {
                "start": row["start"],
                "stop": row["stop"],
                "width": row["width"],
                "constant": f(row["constant"]),
                "fractionOfFull": f(row["constant"] / data["gamma2Full"]) if data["gamma2Full"] else 0.0,
                "lambdaStart": f(row["lambdaStart"]),
                "lambdaStop": f(row["lambdaStop"]),
            }
            for row in data["topWindows"]
        ],
        "topSingles": [
            {
                "mode": row["mode"],
                "constant": f(row["constant"]),
                "fractionOfFull": f(row["constant"] / data["gamma2Full"]) if data["gamma2Full"] else 0.0,
                "lambda": f(row["lambda"]),
            }
            for row in data["topSingles"]
        ],
        "tails": [
            {
                "start": row["start"],
                "constant": f(row["constant"]),
                "fractionOfFull": f(row["constant"] / data["gamma2Full"]) if data["gamma2Full"] else 0.0,
                "lambdaStart": f(row["lambdaStart"]),
            }
            for row in data["tails"]
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
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--json-out", default="windowed_hardy_green_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print(
        f"Windowed Hardy/Green scan model={args.model} kind={args.kind} "
        f"omega={args.omega} L={args.L}",
        flush=True,
    )
    print(
        "  basis cons gamma2      worst window        frac   best single       frac",
        flush=True,
    )
    for basis in parse_ints(args.basis):
        case = make_case(args, basis)
        row = serial_case(case)
        rows.append(row)
        top_window = row["topWindows"][0]
        top_single = row["topSingles"][0]
        wlabel = f"[{top_window['start']},{top_window['stop']})"
        slabel = str(top_single["mode"])
        print(
            f"  {basis:5d} {case.constraints:4d} {fmt(row['gamma2Full']):>10} "
            f"{wlabel:>14} {top_window['fractionOfFull']:7.3f} "
            f"{slabel:>13} {top_single['fractionOfFull']:7.3f}",
            flush=True,
        )

    out = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
