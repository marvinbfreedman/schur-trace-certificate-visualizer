#!/usr/bin/env python3
"""High-precision row-null quadratic check for local normalized nodes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_bordered_minors_mp import beta, p0  # noqa: E402


def rows(c, lam):
    e = [1 / x for x in lam]
    d = [beta(c, x) / (x - c) for x in lam]
    return e, d


def null_basis(e, d):
    # Use columns 0,1 as pivots.  For each free index k>=2, solve for a0,a1
    # with a_k=1 and row constraints zero.
    det = e[0] * d[1] - e[1] * d[0]
    basis = []
    n = len(e)
    for k in range(2, n):
        rhs0 = -e[k]
        rhs1 = -d[k]
        a0 = (rhs0 * d[1] - e[1] * rhs1) / det
        a1 = (e[0] * rhs1 - rhs0 * d[0]) / det
        col = [mp.mpf("0")] * n
        col[0] = a0
        col[1] = a1
        col[k] = mp.mpf("1")
        basis.append(col)
    return basis


def khat(c, x, y):
    return p0(c, x, y) / ((x - c) * (y - c))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h = mp.mpf(args.h)
    s = [s0 + i * h for i in range(args.n)]
    lam = [c * mp.e**x for x in s]
    e, d = rows(c, lam)
    basis = null_basis(e, d)
    m = len(basis)
    q = mp.matrix(m)
    for a in range(m):
        for b in range(m):
            total = mp.mpf("0")
            for i in range(args.n):
                for j in range(args.n):
                    total += basis[a][i] * basis[b][j] * khat(c, lam[i], lam[j])
            q[a, b] = total
    eigvals, eigvecs = mp.eigsy(q, eigvals_only=False)
    detq = mp.det(q)
    print(
        f"endpoint row-null local mp n={args.n} s0={s0} h={h} dps={args.dps}"
    )
    print(f"  det row-null Q = {mp.nstr(detq, 20)}")
    print("  eigvals:")
    for val in eigvals:
        print(f"    {mp.nstr(val, 25)}")
    min_index = min(range(len(eigvals)), key=lambda idx: eigvals[idx])
    if eigvals[min_index] < 0:
        coeff = [mp.mpf("0")] * args.n
        for a in range(m):
            for i in range(args.n):
                coeff[i] += basis[a][i] * eigvecs[a, min_index]
        scale = max(abs(x) for x in coeff)
        coeff = [x / scale for x in coeff]
        print("  negative witness alpha scaled:")
        for i, val in enumerate(coeff):
            print(f"    i={i} s={mp.nstr(s[i], 8)} alpha={mp.nstr(val, 18)}")
    print("  row residual max:")
    max_res = mp.mpf("0")
    for col in basis:
        max_res = max(max_res, abs(sum(ei * ci for ei, ci in zip(e, col))))
        max_res = max(max_res, abs(sum(di * ci for di, ci in zip(d, col))))
    print(f"    {mp.nstr(max_res, 12)}")


if __name__ == "__main__":
    main()
