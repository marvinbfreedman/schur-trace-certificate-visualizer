#!/usr/bin/env python3
import argparse
import itertools
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_tp_check.py requires numpy; run with python") from exc


def h_value(c, u):
    eu = math.exp(u)
    return math.exp(1.25 * u) * (c * eu - 1.5) * math.exp(-c * eu)


def p_values(c, grid, r_end, r_order):
    r_pts, r_wts = quadrature(r_end, r_order)
    out = {}
    for t in grid:
        for s in grid:
            total = 0.0
            for r, w in zip(r_pts, r_wts):
                total += (
                    float(w)
                    * (float(r) + 0.5 * (float(t) + float(s)))
                    * h_value(c, float(t + r))
                    * h_value(c, float(s + r))
                )
            out[(float(t), float(s))] = total
    return out


def scan_h(c, grid, max_order, sample_limit, tol):
    print("  sign-regular minors for h(t+s):")
    for n in range(1, max_order + 1):
        expected = -1 if (n * (n - 1) // 2) % 2 else 1
        worst = float("inf")
        wrong = 0
        total = 0
        for rows in itertools.combinations(grid, n):
            if total >= sample_limit:
                break
            for cols in itertools.combinations(grid, n):
                mat = np.array(
                    [[h_value(c, float(t + s)) for s in cols] for t in rows],
                    dtype=float,
                )
                value = expected * float(np.linalg.det(mat))
                worst = min(worst, value)
                wrong += value < -tol
                total += 1
                if total >= sample_limit:
                    break
        print(
            f"    order={n}: expected_sign={expected:+d} "
            f"min_signed={worst:.12e} wrong={wrong}/{total}"
        )


def scan_p(c, grid, max_order, sample_limit, tol, r_end, r_order):
    values = p_values(c, grid, r_end, r_order)
    print("  total-positive minors for P(t,s):")
    for n in range(1, max_order + 1):
        worst = float("inf")
        wrong = 0
        total = 0
        for rows in itertools.combinations(grid, n):
            if total >= sample_limit:
                break
            for cols in itertools.combinations(grid, n):
                mat = np.array(
                    [[values[(float(t), float(s))] for s in cols] for t in rows],
                    dtype=float,
                )
                value = float(np.linalg.det(mat))
                worst = min(worst, value)
                wrong += value < -tol
                total += 1
                if total >= sample_limit:
                    break
        print(f"    order={n}: min={worst:.12e} wrong={wrong}/{total}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument("--grid", type=int, default=14)
    parser.add_argument("--max-order", type=int, default=5)
    parser.add_argument("--sample-limit", type=int, default=2000)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--r-order", type=int, default=260)
    parser.add_argument("--tol", type=float, default=1e-14)
    args = parser.parse_args()

    grid = np.linspace(0.0, args.T, args.grid)
    print(
        f"endpoint TP check grid=[0,{args.T:g}] grid_n={args.grid} "
        f"max_order={args.max_order}"
    )
    scan_h(args.c, grid, args.max_order, args.sample_limit, args.tol)
    scan_p(
        args.c,
        grid,
        args.max_order,
        args.sample_limit,
        args.tol,
        args.rmax,
        args.r_order,
    )


if __name__ == "__main__":
    main()
