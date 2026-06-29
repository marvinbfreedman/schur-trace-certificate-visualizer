#!/usr/bin/env python3
import argparse
import math

from klm_test import PI


def ranged_phi_derivatives_abs(x: float, nmin: int, nmax: int):
    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
    total2 = 0.0
    for n in range(nmin, nmax + 1):
        n2 = n * n
        c = PI * n2
        for coeff, lam in (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        ):
            exp_value = math.exp(lam * x - c * e2)
            log_derivative = lam - 2.0 * c * e2
            second_factor = log_derivative * log_derivative - 4.0 * c * e2
            term0 = coeff * exp_value
            total0 += term0
            total1 += term0 * log_derivative
            total2 += term0 * second_factor
    return total0, total1, total2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nmin", type=int, default=4)
    parser.add_argument("--nmax", type=int, default=40)
    parser.add_argument("--xmax", type=float, default=8.0)
    parser.add_argument("--n-grid", type=int, default=2000)
    args = parser.parse_args()

    max0 = (0.0, None)
    max1 = (0.0, None)
    max2 = (0.0, None)
    for i in range(args.n_grid):
        x = args.xmax * i / (args.n_grid - 1)
        f0, f1, f2 = ranged_phi_derivatives_abs(x, args.nmin, args.nmax)
        if abs(f0) > max0[0]:
            max0 = (abs(f0), x)
        if abs(f1) > max1[0]:
            max1 = (abs(f1), x)
        if abs(f2) > max2[0]:
            max2 = (abs(f2), x)

    print(f"tail n=[{args.nmin},{args.nmax}] x=[0,{args.xmax:g}] grid={args.n_grid}")
    print(f"  sup |tail|  ~= {max0[0]:.12e} at x={max0[1]:.6g}")
    print(f"  sup |tail'| ~= {max1[0]:.12e} at x={max1[1]:.6g}")
    print(f"  sup |tail''|~= {max2[0]:.12e} at x={max2[1]:.6g}")


if __name__ == "__main__":
    main()
