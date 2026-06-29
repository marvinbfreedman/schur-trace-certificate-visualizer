#!/usr/bin/env python3
"""Bordered determinant tests for normalized conditional positivity.

For rows e,d and kernel Phat, conditional positive definiteness on
ker(e,d) implies a sign condition for the bordered matrix

  [ 0 0 | e(s_j) ]
  [ 0 0 | d(s_j) ]
  [ e(s_i) d(s_i) | Phat(s_i,s_j) ].

With two constraints the expected nondegenerate determinant sign is positive.
This script tests consecutive-grid and random minors.  Stable positivity would
suggest the next analytic theorem: bordered total positivity of
(e,d,Phat).
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import beta_closed, c_kernel, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_bordered_minors.py requires numpy") from exc


def khat(c: float, s: float, t: float, kernel) -> float:
    lam = c * math.exp(s)
    mu = c * math.exp(t)
    return kernel(c, lam, mu) / ((lam - c) * (mu - c))


def rows(c: float, s: float) -> tuple[float, float]:
    lam = c * math.exp(s)
    e = 1.0 / lam
    d = beta_closed(c, lam) / (lam - c)
    return e, d


def bordered_slogdet(c: float, s_pts, kernel) -> tuple[float, float]:
    n = len(s_pts)
    mat = np.zeros((n + 2, n + 2), dtype=float)
    for j, s in enumerate(s_pts):
        e, d = rows(c, float(s))
        mat[0, j + 2] = mat[j + 2, 0] = e
        mat[1, j + 2] = mat[j + 2, 1] = d
    for i, s in enumerate(s_pts):
        for j, t in enumerate(s_pts[: i + 1]):
            val = khat(c, float(s), float(t), kernel)
            mat[i + 2, j + 2] = mat[j + 2, i + 2] = val
    sign, logabs = np.linalg.slogdet(mat)
    return float(sign), float(logabs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--grid", type=int, default=30)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--trials", type=int, default=300)
    parser.add_argument("--seed", type=int, default=99)
    parser.add_argument("--start", type=float, default=1e-3)
    parser.add_argument("--kernel", choices=("p0", "c"), default="p0")
    args = parser.parse_args()

    kernel = p0_kernel if args.kernel == "p0" else c_kernel
    rng = np.random.default_rng(args.seed)
    grid = np.linspace(args.start, args.T, args.grid)
    print(
        f"endpoint normalized bordered minors kernel={args.kernel} "
        f"s=[{args.start:g},{args.T:g}] "
        f"grid={args.grid} trials={args.trials}"
    )
    for n in range(2, args.max_n + 1):
        signs = []
        logs = []
        for left in range(0, args.grid - n + 1):
            sign, logabs = bordered_slogdet(args.c, grid[left : left + n], kernel)
            signs.append(sign)
            logs.append(logabs)
        for _ in range(args.trials):
            pts = np.sort(rng.uniform(args.start, args.T, size=n))
            sign, logabs = bordered_slogdet(args.c, pts, kernel)
            signs.append(sign)
            logs.append(logabs)
        signs = np.array(signs)
        logs = np.array(logs)
        bad = int((signs < 0).sum())
        zero = int((signs == 0).sum())
        print(
            f"  n={n:2d}: sign_bad={bad:4d} zero={zero:4d} "
            f"logabs_min={logs.min(): .6e} logabs_max={logs.max(): .6e}"
        )


if __name__ == "__main__":
    main()
