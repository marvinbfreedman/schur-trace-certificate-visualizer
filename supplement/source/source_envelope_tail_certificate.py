#!/usr/bin/env python3
r"""Source-envelope certificate for high-frequency tail decay.

Let A e_j=lambda_j e_j and b_j(u)=<B u,e_j>.  Then

    C_tail(M)
      = || sum_{j>=M} lambda_j^{-1} b_j^* b_j ||
      <= sum_{j>=M} lambda_j^{-1} ||b_j||^2.

The scalar envelope

    s_j = lambda_j^{-1} ||b_j||^2 / C_full

therefore gives a direct sufficient condition for tail decay:

    sum_{j>=M} s_j -> 0.

This is sharper than generic spectral moments because it keeps the actual
source projections.
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

from high_frequency_tail_refinement import constraints_for  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, max_eig_or_zero, trace_matrix  # noqa: E402


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def f(x):
    return float(x)


def fmt(x, digits=8):
    return mp.nstr(x, digits)


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


def compute_case(base, basis):
    constraints = constraints_for(base, basis)
    args = make_args(base, basis, constraints)
    K, polys = gram_matrix(args, mp.mpf(args.omega), mp.mpf(args.L))
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
    beig = avecs.T * B
    gamma2 = window_constant(beig, avals, 0, len(avals))

    envelope = []
    for j, lam in enumerate(avals):
        row2 = row_norm2(beig, j)
        scalar = row2 / lam
        envelope.append(
            {
                "mode": j,
                "lambda": lam,
                "rowNorm2": row2,
                "scalar": scalar,
                "scalarFrac": scalar / gamma2 if gamma2 else mp.mpf("0"),
                "operatorSingleFrac": window_constant(beig, avals, j, j + 1) / gamma2 if gamma2 else mp.mpf("0"),
            }
        )

    tail_rows = []
    for start in range(len(avals) + 1):
        op_tail = window_constant(beig, avals, start, len(avals))
        scalar_tail = mp.fsum(item["scalar"] for item in envelope[start:])
        tail_rows.append(
            {
                "start": start,
                "operatorTail": op_tail,
                "operatorFrac": op_tail / gamma2 if gamma2 else mp.mpf("0"),
                "scalarTail": scalar_tail,
                "scalarFrac": scalar_tail / gamma2 if gamma2 else mp.mpf("0"),
                "slack": op_tail / scalar_tail if scalar_tail else mp.mpf("0"),
            }
        )

    return {
        "basis": basis,
        "constraints": constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "gamma2": gamma2,
        "envelope": envelope,
        "tails": tail_rows,
    }


def serial(row):
    return {
        "basis": row["basis"],
        "constraints": row["constraints"],
        "rank": row["rank"],
        "nullity": row["nullity"],
        "positiveModes": row["positiveModes"],
        "gamma2": f(row["gamma2"]),
        "envelope": [
            {
                "mode": item["mode"],
                "lambda": f(item["lambda"]),
                "rowNorm2": f(item["rowNorm2"]),
                "scalar": f(item["scalar"]),
                "scalarFrac": f(item["scalarFrac"]),
                "operatorSingleFrac": f(item["operatorSingleFrac"]),
            }
            for item in row["envelope"]
        ],
        "tails": [
            {
                "start": item["start"],
                "operatorTail": f(item["operatorTail"]),
                "operatorFrac": f(item["operatorFrac"]),
                "scalarTail": f(item["scalarTail"]),
                "scalarFrac": f(item["scalarFrac"]),
                "slack": f(item["slack"]),
            }
            for item in row["tails"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="12 14 16 18 20")
    parser.add_argument("--constraint-rule", choices=["ratio", "target-nullity", "offset"], default="ratio")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--target-nullity", type=int, default=6)
    parser.add_argument("--constraint-offset", type=int, default=6)
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
    parser.add_argument("--json-out", default="source_envelope_tail_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print(
        f"Source envelope tail certificate model={args.model} kind={args.kind}",
        flush=True,
    )
    print("  basis cons null gamma2      tail4_op tail4_scalar tail6_op tail6_scalar", flush=True)
    for basis in parse_ints(args.basis):
        row = compute_case(args, basis)
        if row is None:
            continue
        rows.append(serial(row))
        tails = row["tails"]
        def frac_at(idx, key):
            if idx >= len(tails):
                return mp.mpf("0")
            return tails[idx][key]
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['nullity']:4d} "
            f"{fmt(row['gamma2']):>10} "
            f"{fmt(frac_at(4, 'operatorFrac')):>9} {fmt(frac_at(4, 'scalarFrac')):>12} "
            f"{fmt(frac_at(6, 'operatorFrac')):>9} {fmt(frac_at(6, 'scalarFrac')):>12}",
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
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
