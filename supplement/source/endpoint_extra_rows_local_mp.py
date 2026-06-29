#!/usr/bin/env python3
"""High-precision local row-null test with optional extra rows."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_bordered_minors_mp import beta, p0  # noqa: E402


def row_value(name, c, lam):
    if name == "e":
        return 1 / lam
    if name == "d":
        return beta(c, lam) / (lam - c)
    if name == "r":
        return mp.log(lam / c) / (lam - c)
    if name == "one":
        return mp.mpf("1")
    if name == "s":
        return mp.log(lam / c)
    if name == "lambda":
        return lam
    raise ValueError(f"unknown row {name}")


def null_basis(rows):
    m = len(rows)
    n = len(rows[0])
    pivot = mp.matrix([[rows[a][b] for b in range(m)] for a in range(m)])
    basis = []
    for k in range(m, n):
        rhs = mp.matrix([-rows[a][k] for a in range(m)])
        sol = mp.lu_solve(pivot, rhs)
        col = [mp.mpf("0")] * n
        for i in range(m):
            col[i] = sol[i]
        col[k] = mp.mpf("1")
        basis.append(col)
    return basis


def khat(c, x, y):
    return p0(c, x, y) / ((x - c) * (y - c))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", default="e,d,r")
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=120)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    row_names = [x.strip() for x in args.rows.split(",") if x.strip()]
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h = mp.mpf(args.h)
    s = [s0 + i * h for i in range(args.n)]
    lam = [c * mp.e**x for x in s]
    rows = [[row_value(name, c, x) for x in lam] for name in row_names]
    basis = null_basis(rows)
    m = len(basis)
    q = mp.matrix(m)
    for a in range(m):
        for b in range(m):
            total = mp.mpf("0")
            for i in range(args.n):
                for j in range(args.n):
                    total += basis[a][i] * basis[b][j] * khat(c, lam[i], lam[j])
            q[a, b] = total
    eigvals = mp.eigsy(q, eigvals_only=True) if m else []
    print(
        f"endpoint extra rows local mp rows={','.join(row_names)} "
        f"n={args.n} s0={s0} h={h} dps={args.dps}"
    )
    if m:
        print(f"  det Q = {mp.nstr(mp.det(q), 20)}")
        print("  eigvals:")
        for val in eigvals:
            print(f"    {mp.nstr(val, 25)}")
    else:
        print("  no nullspace")


if __name__ == "__main__":
    main()
