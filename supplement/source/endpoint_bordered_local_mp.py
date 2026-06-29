#!/usr/bin/env python3
"""High-precision local bordered determinant checks.

For a totally positive/bordered-TP kernel, determinants should remain positive
as nodes coalesce.  This script evaluates the bordered determinant on

  s_i=s0+i*h, i=0,...,n-1,

with mpmath arithmetic.  Positivity here is evidence for a confluent
Chebyshev/Wronskian theorem behind the bordered determinant route.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_bordered_minors_mp import bordered_det  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--h", type=str, default="0.02")
    parser.add_argument("--s0-min", type=str, default="0.001")
    parser.add_argument("--s0-max", type=str, default="8")
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--dps", type=int, default=90)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    h = mp.mpf(args.h)
    s0_min = mp.mpf(args.s0_min)
    s0_max = mp.mpf(args.s0_max)
    bad = 0
    smallest = None
    smallest_s0 = None
    print(
        f"endpoint bordered local mp n={args.n} h={h} "
        f"s0=[{s0_min},{s0_max}] samples={args.samples} dps={args.dps}"
    )
    for k in range(args.samples):
        if args.samples == 1:
            s0 = s0_min
        else:
            s0 = s0_min + (s0_max - s0_min) * k / (args.samples - 1)
        pts = [s0 + i * h for i in range(args.n)]
        det = bordered_det(c, pts)
        if det < 0:
            bad += 1
        if smallest is None or det < smallest:
            smallest = det
            smallest_s0 = s0
        print(
            f"  s0={mp.nstr(s0, 8)} det={mp.nstr(det, 10)} "
            f"log10={mp.nstr(mp.log10(abs(det)), 8)}"
        )
    print(f"  bad signs = {bad}")
    print(
        f"  smallest at s0={mp.nstr(smallest_s0, 10)} "
        f"det={mp.nstr(smallest, 20)}"
    )


if __name__ == "__main__":
    main()
