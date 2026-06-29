#!/usr/bin/env python3
"""Finite-difference complete-monotonicity scan for d(s)=beta/(lambda-c).

If d(s) is a Laplace transform of a positive measure, then uniform forward
differences satisfy

  (-1)^k Delta^k d >= 0

for every k.  This is only a diagnostic, but it quickly rules out the
complete-monotone proof route if a robust sign failure appears.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import beta_closed

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_d_cm_scan.py requires numpy") from exc


def d_value(c: float, s: float) -> float:
    lam = c * math.exp(s)
    return beta_closed(c, lam) / (lam - c)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=12.0)
    parser.add_argument("--points", type=int, default=400)
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--start", type=float, default=1e-4)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    s = np.linspace(args.start, args.T, args.points)
    diff = np.array([d_value(args.c, float(x)) for x in s])
    print(
        f"endpoint d complete-monotone finite differences "
        f"s=[{args.start:g},{args.T:g}] points={args.points}"
    )
    for k in range(args.max_order + 1):
        signed = ((-1.0) ** k) * diff
        print(
            f"  order {k:2d}: min signed={signed.min(): .12e} "
            f"max signed={signed.max(): .12e} "
            f"bad={(signed < -args.tol).sum():4d}"
        )
        diff = np.diff(diff)


if __name__ == "__main__":
    main()
