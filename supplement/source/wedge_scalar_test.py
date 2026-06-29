#!/usr/bin/env python3
import argparse
import math

from klm_test import phi, simpson_grid


def log_phi(x: float) -> float:
    value = phi(x)
    if value <= 0.0:
        return float("-inf")
    return math.log(value)


def score(x: float, step: float) -> float:
    plus = log_phi(x + step)
    minus = log_phi(x - step)
    if not math.isfinite(plus) or not math.isfinite(minus):
        return float("nan")
    return -(plus - minus) / (2.0 * step)


def weight(z: float, omega: float) -> float:
    return 0.25 * z * math.cosh(omega * z)


def reflected_source(c: float, v: float, omega: float) -> float:
    return weight(c - 2.0 * v, omega) * phi(c - v) * phi(v)


def anti_diagonal_derivative(c: float, v: float, omega: float, rs, ws, step: float) -> float:
    total = 0.0
    for r, w in zip(rs, ws):
        left = c - v + r
        right = v + r
        left_phi = phi(left)
        right_phi = phi(right)
        if left_phi <= 0.0 or right_phi <= 0.0:
            continue
        score_gap = score(left, step) - score(right, step)
        if not math.isfinite(score_gap):
            continue
        total += (
            w
            * weight(c + 2.0 * r, omega)
            * left_phi
            * right_phi
            * score_gap
        )
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--cmax", type=float, default=8.0)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=900)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    rs, ws = simpson_grid(args.rmax, args.intervals)
    min_gap = float("inf")
    min_at = None
    max_ratio = -float("inf")
    max_ratio_at = None

    for i in range(1, args.n):
        c = args.cmax * i / (args.n - 1)
        for j in range(i // 2 + 1):
            v = 0.5 * c * j / max(1, i // 2)
            lhs = anti_diagonal_derivative(c, v, args.omega, rs, ws, args.h)
            rhs = reflected_source(c, v, args.omega)
            gap = lhs - rhs
            if gap < min_gap:
                min_gap = gap
                min_at = (c, v, lhs, rhs)
            if lhs > 0:
                ratio = rhs / lhs
                if ratio > max_ratio:
                    max_ratio = ratio
                    max_ratio_at = (c, v)

    print(f"omega={args.omega:g} c in [0,{args.cmax:g}] n={args.n}")
    print(
        "  min anti-diagonal gap="
        f"{min_gap:.12e} at c={min_at[0]:.6g}, v={min_at[1]:.6g} "
        f"lhs={min_at[2]:.12e} rhs={min_at[3]:.12e}"
    )
    print(f"  max rhs/lhs={max_ratio:.12e} at c={max_ratio_at[0]:.6g}, v={max_ratio_at[1]:.6g}")


if __name__ == "__main__":
    main()
