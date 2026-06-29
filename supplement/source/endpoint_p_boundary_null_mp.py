#!/usr/bin/env python3
"""High-precision P-kernel restricted to the endpoint boundary-null row."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_p_kernel_closed_mp import p_entry  # noqa: E402


def boundary_row(c, t):
    return mp.e ** (mp.mpf("1.5") * t - c * mp.e**t)


def null_basis(row):
    # Pivot on first entry; free coordinates 1..n-1.
    basis = []
    n = len(row)
    for k in range(1, n):
        col = [mp.mpf("0")] * n
        col[0] = -row[k] / row[0]
        col[k] = mp.mpf("1")
        basis.append(col)
    return basis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=120)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    step = mp.mpf(args.h)
    pts = [s0 + i * step for i in range(args.n)]
    pmat = mp.matrix(args.n)
    for i in range(args.n):
        for j in range(i + 1):
            val = p_entry(c, pts[i], pts[j])
            pmat[i, j] = pmat[j, i] = val
    full_vals = mp.eigsy(pmat, eigvals_only=True)
    row = [boundary_row(c, t) for t in pts]
    basis = null_basis(row)
    m = len(basis)
    q = mp.matrix(m)
    for a in range(m):
        for b in range(m):
            total = mp.mpf("0")
            for i in range(args.n):
                for j in range(args.n):
                    total += basis[a][i] * basis[b][j] * pmat[i, j]
            q[a, b] = total
    null_vals = mp.eigsy(q, eigvals_only=True)
    print(f"endpoint P boundary-null mp n={args.n} s0={s0} h={step} dps={args.dps}")
    print("  full eig low:")
    for val in full_vals[: min(args.n, 6)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  boundary-null eig low:")
    for val in null_vals[: min(m, 8)]:
        print(f"    {mp.nstr(val, 25)}")


if __name__ == "__main__":
    main()
