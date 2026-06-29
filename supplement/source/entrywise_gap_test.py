#!/usr/bin/env python3
import argparse

from kernel_parity_test import kernel_pm, kernel_pp
from klm_test import simpson_grid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmax", type=float, default=8.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=900)
    args = parser.parse_args()

    xs = [args.xmax * i / (args.n - 1) for i in range(args.n)]
    rs, ws = simpson_grid(args.rmax, args.intervals)
    min_odd = float("inf")
    min_even = float("inf")
    min_odd_at = None
    min_even_at = None
    max_ratio = -float("inf")
    max_ratio_at = None

    for x in xs:
        for y in xs:
            a = kernel_pp(x, y, args.omega, rs, ws)
            b = kernel_pm(x, y, args.omega, rs, ws)
            odd = a - b
            even = a + b
            if odd < min_odd:
                min_odd = odd
                min_odd_at = (x, y)
            if even < min_even:
                min_even = even
                min_even_at = (x, y)
            if a > 0:
                ratio = abs(b) / a
                if ratio > max_ratio:
                    max_ratio = ratio
                    max_ratio_at = (x, y)

    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g}")
    print(f"  min entry A-B={min_odd:.12e} at x={min_odd_at[0]:.6g}, y={min_odd_at[1]:.6g}")
    print(f"  min entry A+B={min_even:.12e} at x={min_even_at[0]:.6g}, y={min_even_at[1]:.6g}")
    print(f"  max entry |B|/A={max_ratio:.12e} at x={max_ratio_at[0]:.6g}, y={max_ratio_at[1]:.6g}")


if __name__ == "__main__":
    main()
