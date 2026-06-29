#!/usr/bin/env python3
"""High-precision bordered determinant checks using local mpmath.

This is a targeted follow-up for endpoint_normalized_bordered_minors.py.  The
double-precision determinants become as small as exp(-200), so isolated signs
past order 9 are not trustworthy.  This script recomputes selected orders with
mpmath arithmetic.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402


def j_values(a):
    root = mp.sqrt(a)
    j_half = 1 / a + mp.sqrt(mp.pi) * mp.e**a * mp.erfc(root) / (2 * a ** mp.mpf("1.5"))
    j_3_2 = 1 / a + mp.mpf("1.5") * j_half / a
    j_5_2 = 1 / a + mp.mpf("2.5") * j_3_2 / a
    j_7_2 = 1 / a + mp.mpf("3.5") * j_5_2 / a
    return j_3_2, j_5_2, j_7_2


def g_gram(lam, mu):
    a = lam + mu
    j3, j5, j7 = j_values(a)
    return j3 - a * j5 + lam * mu * j7


def e_kernel(c, lam, mu):
    return g_gram(lam, mu) - g_gram(lam, c) - g_gram(c, mu) + g_gram(c, c)


def beta(c, lam):
    return g_gram(lam, c) - g_gram(c, c)


def p0(c, lam, mu):
    return mp.mpf("0.5") * mp.log(lam * mu / (c * c)) * e_kernel(c, lam, mu)


def bordered_det(c, s_pts):
    n = len(s_pts)
    mat = mp.matrix(n + 2)
    for j, s0 in enumerate(s_pts):
        s = mp.mpf(str(float(s0)))
        lam = c * mp.e**s
        e = 1 / lam
        d = beta(c, lam) / (lam - c)
        mat[0, j + 2] = mat[j + 2, 0] = e
        mat[1, j + 2] = mat[j + 2, 1] = d
    for i, s0 in enumerate(s_pts):
        s = mp.mpf(str(float(s0)))
        lam = c * mp.e**s
        for j, t0 in enumerate(s_pts[: i + 1]):
            t = mp.mpf(str(float(t0)))
            mu = c * mp.e**t
            val = p0(c, lam, mu) / ((lam - c) * (mu - c))
            mat[i + 2, j + 2] = mat[j + 2, i + 2] = val
    return mp.det(mat)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--start", type=float, default=1e-3)
    parser.add_argument("--grid", type=int, default=34)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--trials", type=int, default=30)
    parser.add_argument("--seed", type=int, default=99)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    rng = np.random.default_rng(args.seed)
    grid = np.linspace(args.start, args.T, args.grid)
    bad = 0
    smallest = None
    smallest_pts = None

    samples = [grid[left : left + args.n] for left in range(0, args.grid - args.n + 1)]
    samples.extend(
        np.sort(rng.uniform(args.start, args.T, size=args.n))
        for _ in range(args.trials)
    )
    for pts in samples:
        det = bordered_det(c, pts)
        if det < 0:
            bad += 1
        if smallest is None or det < smallest:
            smallest = det
            smallest_pts = pts

    print(
        f"endpoint bordered mp n={args.n} T={args.T:g} samples={len(samples)} "
        f"dps={args.dps}"
    )
    print(f"  bad signs = {bad}")
    print(f"  smallest det = {mp.nstr(smallest, 20)}")
    print(f"  log10(|smallest|) = {mp.nstr(mp.log10(abs(smallest)), 12)}")
    print(
        "  smallest pts = "
        + " ".join(f"{float(x):.6g}" for x in smallest_pts[: min(len(smallest_pts), 12)])
    )


if __name__ == "__main__":
    main()
