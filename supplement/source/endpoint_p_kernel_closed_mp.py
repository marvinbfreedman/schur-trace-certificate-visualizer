#!/usr/bin/env python3
"""High-precision closed-form local tests for the endpoint P kernel.

This uses

  int_1^inf x^p exp(-A x) dx = A^(-p-1) Gamma(p+1,A)

and differentiates in p to obtain the log(x) moment.  It avoids numerical
quadrature in the highly confluent determinant tests for P(t,s).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


def ipow(p, a):
    return a ** (-(p + 1)) * mp.gammainc(p + 1, a, mp.inf)


def ilog(p, a):
    return mp.diff(lambda pp: ipow(pp, a), p)


def p_entry(c, t, s):
    et = mp.e**t
    es = mp.e**s
    a = c * (et + es)
    pref = mp.e ** (mp.mpf("1.25") * (t + s))
    const = (t + s) / 2
    coeffs = [
        (mp.mpf("3.5"), c * c * et * es),
        (mp.mpf("2.5"), -mp.mpf("1.5") * c * (et + es)),
        (mp.mpf("1.5"), mp.mpf("2.25")),
    ]
    total = mp.mpf("0")
    for p, coeff in coeffs:
        total += coeff * (ilog(p, a) + const * ipow(p, a))
    return pref * total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    step = mp.mpf(args.h)
    pts = [s0 + i * step for i in range(args.n)]
    mat = mp.matrix(args.n)
    for i in range(args.n):
        for j in range(i + 1):
            val = p_entry(c, pts[i], pts[j])
            mat[i, j] = mat[j, i] = val
    eigvals = mp.eigsy(mat, eigvals_only=True)
    det = mp.det(mat)
    print(f"endpoint P closed mp n={args.n} s0={s0} h={step} dps={args.dps}")
    print(f"  det = {mp.nstr(det, 25)}")
    print("  eig low:")
    for val in eigvals[: min(args.n, 10)]:
        print(f"    {mp.nstr(val, 25)}")


if __name__ == "__main__":
    main()
