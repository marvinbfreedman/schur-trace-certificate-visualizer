#!/usr/bin/env python3
"""Verify the exact bordered determinant invariance Chat vs Phat.

After endpoint normalization,

  Chat-Phat = 1/2[d(lambda) r(mu)+r(lambda)d(mu)],
  r(lambda)=log(lambda/c)/(lambda-c).

Since d is one of the two bordered rows, adding this symmetric row-containing
rank term leaves the row-null quadratic form and bordered determinants
unchanged.  This script checks the algebra numerically on random nodes.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import beta_closed, c_kernel, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_chat_bordered_invariance.py requires numpy") from exc


def normed_rows(c, lam):
    d = np.array([beta_closed(c, float(x)) / (float(x) - c) for x in lam])
    r = np.array([math.log(float(x) / c) / (float(x) - c) for x in lam])
    return d, r


def bordered(c, lam, kernel):
    n = len(lam)
    mat = np.zeros((n + 2, n + 2), dtype=float)
    d, _ = normed_rows(c, lam)
    e = 1.0 / lam
    mat[0, 2:] = mat[2:, 0] = e
    mat[1, 2:] = mat[2:, 1] = d
    for i, x in enumerate(lam):
        for j, y in enumerate(lam[: i + 1]):
            val = kernel(c, float(x), float(y)) / ((float(x) - c) * (float(y) - c))
            mat[i + 2, j + 2] = mat[j + 2, i + 2] = val
    return mat


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    s = np.sort(rng.uniform(1e-3, args.T, size=args.n))
    lam = args.c * np.exp(s)
    d, r = normed_rows(args.c, lam)
    ph = bordered(args.c, lam, p0_kernel)
    ch = bordered(args.c, lam, c_kernel)
    expected = np.zeros_like(ph)
    expected[2:, 2:] = 0.5 * (np.outer(d, r) + np.outer(r, d))
    print(f"endpoint Chat bordered invariance n={args.n} T={args.T:g}")
    print(f"  ||Chat-Phat-rankterm|| = {np.linalg.norm(ch-ph-expected):.12e}")
    print(f"  det Phat bordered      = {np.linalg.det(ph):.12e}")
    print(f"  det Chat bordered      = {np.linalg.det(ch):.12e}")
    print(f"  determinant defect     = {np.linalg.det(ch)-np.linalg.det(ph):.12e}")


if __name__ == "__main__":
    main()
