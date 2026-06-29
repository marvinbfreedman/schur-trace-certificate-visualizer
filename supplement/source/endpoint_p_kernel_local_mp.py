#!/usr/bin/env python3
"""High-precision local tests for the endpoint Green kernel P(t,s).

P(t,s)=int_0^inf (r+(t+s)/2) h(t+r)h(s+r) dr

is the Hankel/Green kernel equivalent to the A-branch endpoint range Hardy
inequality.  This script tests close-node principal matrices, where determinant
and finite-index shortcuts are most likely to fail.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


def h_value(c, u):
    z = c * mp.e**u
    return mp.e ** (mp.mpf("1.25") * u) * (z - mp.mpf("1.5")) * mp.e ** (-z)


def p_entry(c, t, s, rmax):
    def integrand(r):
        return (r + (t + s) / 2) * h_value(c, t + r) * h_value(c, s + r)

    points = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    x = mp.mpf("0.4")
    while x < rmax:
        points.append(x)
        x *= 2
    points.append(rmax)
    return mp.quad(integrand, points)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--rmax", type=str, default="12")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    step = mp.mpf(args.h)
    rmax = mp.mpf(args.rmax)
    pts = [s0 + i * step for i in range(args.n)]
    mat = mp.matrix(args.n)
    for i in range(args.n):
        for j in range(i + 1):
            val = p_entry(c, pts[i], pts[j], rmax)
            mat[i, j] = mat[j, i] = val
    eigvals = mp.eigsy(mat, eigvals_only=True)
    det = mp.det(mat)
    print(f"endpoint P local mp n={args.n} s0={s0} h={step} dps={args.dps}")
    print(f"  det = {mp.nstr(det, 25)}")
    print("  eig low:")
    for val in eigvals[: min(args.n, 10)]:
        print(f"    {mp.nstr(val, 25)}")


if __name__ == "__main__":
    main()
