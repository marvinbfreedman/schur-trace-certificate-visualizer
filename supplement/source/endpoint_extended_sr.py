#!/usr/bin/env python3
import argparse
import itertools
import math

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_extended_sr.py requires numpy; run with python") from exc


def h_value(c, u):
    eu = math.exp(u)
    return math.exp(1.25 * u) * (c * eu - 1.5) * math.exp(-c * eu)


def extended_value(c, t, r, kind):
    u = t + r
    value = h_value(c, u)
    if kind == 0:
        return value
    if kind == 1:
        return u * value
    raise ValueError(kind)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument("--grid", type=int, default=12)
    parser.add_argument("--max-order", type=int, default=5)
    parser.add_argument("--sample-limit", type=int, default=1200)
    parser.add_argument("--tol", type=float, default=1e-14)
    args = parser.parse_args()

    grid = np.linspace(0.0, args.T, args.grid)
    print(
        f"endpoint extended sign-regularity grid=[0,{args.T:g}] "
        f"grid_n={args.grid} max_order={args.max_order}"
    )
    for n in range(1, args.max_order + 1):
        expected = -1 if (n * (n - 1) // 2) % 2 else 1
        print(f"  order={n} expected_sign={expected:+d}")
        for kinds in itertools.product((0, 1), repeat=n):
            worst = float("inf")
            wrong = 0
            total = 0
            for rows in itertools.combinations(grid, n):
                if total >= args.sample_limit:
                    break
                for cols in itertools.combinations(grid, n):
                    mat = np.array(
                        [
                            [
                                extended_value(args.c, float(t), float(r), kinds[j])
                                for j, r in enumerate(cols)
                            ]
                            for t in rows
                        ],
                        dtype=float,
                    )
                    value = expected * float(np.linalg.det(mat))
                    worst = min(worst, value)
                    wrong += value < -args.tol
                    total += 1
                    if total >= args.sample_limit:
                        break
            label = "".join(str(k) for k in kinds)
            print(
                f"    kinds={label:<{n}} min_signed={worst:.12e} "
                f"wrong={wrong}/{total}"
            )


if __name__ == "__main__":
    main()
