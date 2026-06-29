#!/usr/bin/env python3
"""Random row-null tests for the normalized range theorem.

Build random coefficient vectors alpha, project them onto the nullspace of the
two normalized rows

  e_i=1/lambda_i,  d_i=beta(lambda_i)/(lambda_i-c),

and evaluate Phat0(alpha).  This is a direct sanity check of the conditional
positivity theorem in the source-moment normalization.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import beta_closed, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_row_null_random.py requires numpy") from exc


def project_null(alpha, rows):
    r = np.vstack(rows)
    gram = r @ r.T
    return alpha - r.T @ np.linalg.solve(gram, r @ alpha)


def quadratic(c, lam, alpha):
    total = 0.0
    for i, li in enumerate(lam):
        for j, lj in enumerate(lam):
            total += (
                alpha[i]
                * alpha[j]
                * p0_kernel(c, float(li), float(lj))
                / ((li - c) * (lj - c))
            )
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--n", type=int, default=30)
    parser.add_argument("--trials", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    s = np.linspace(args.T / args.n, args.T, args.n)
    lam = args.c * np.exp(s)
    erow = 1.0 / lam
    drow = np.array(
        [beta_closed(args.c, float(x)) / (float(x) - args.c) for x in lam]
    )
    worst = float("inf")
    worst_rows = None
    for _ in range(args.trials):
        alpha = rng.normal(size=args.n)
        alpha = project_null(alpha, [erow, drow])
        val = quadratic(args.c, lam, alpha)
        if val < worst:
            worst = val
            worst_rows = (float(erow @ alpha), float(drow @ alpha))
    print(
        f"endpoint normalized row-null random T={args.T:g} n={args.n} "
        f"trials={args.trials}"
    )
    print(f"  worst Phat0(alpha) = {worst:.12e}")
    print(f"  worst row residuals e,d = {worst_rows[0]:.3e}, {worst_rows[1]:.3e}")


if __name__ == "__main__":
    main()
