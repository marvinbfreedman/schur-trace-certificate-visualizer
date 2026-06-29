#!/usr/bin/env python3
import argparse
import math

from klm_test import phi


def log_phi(x: float) -> float:
    value = phi(x)
    if value <= 0.0:
        return float("nan")
    return math.log(value)


def score(x: float, h: float) -> float:
    plus = log_phi(x + h)
    minus = log_phi(x - h)
    if not math.isfinite(plus) or not math.isfinite(minus):
        return float("nan")
    return -(plus - minus) / (2.0 * h)


def score_prime(x: float, h: float) -> float:
    center = log_phi(x)
    plus = log_phi(x + h)
    minus = log_phi(x - h)
    if not math.isfinite(center) or not math.isfinite(plus) or not math.isfinite(minus):
        return float("nan")
    return -(plus - 2.0 * center + minus) / (h * h)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xmin", type=float, default=1e-4)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=400)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    xs = [args.xmin + (args.xmax - args.xmin) * i / (args.n - 1) for i in range(args.n)]
    h_values = []
    hp_values = []
    hp_over_h = []
    x_hp_over_h = []
    h_over_x = []
    for x in xs:
        h_value = score(x, args.h)
        hp_value = score_prime(x, args.h)
        h_values.append(h_value)
        hp_values.append(hp_value)
        hp_over_h.append(hp_value / h_value)
        x_hp_over_h.append(x * hp_value / h_value)
        h_over_x.append(h_value / x)

    def min_diff(values):
        best = float("inf")
        best_i = None
        for i in range(len(values) - 1):
            diff = values[i + 1] - values[i]
            if diff < best:
                best = diff
                best_i = i
        return best, best_i

    hp_h_diff, hp_h_i = min_diff(hp_over_h)
    x_hp_h_diff, x_hp_h_i = min_diff(x_hp_over_h)
    h_x_diff, h_x_i = min_diff(h_over_x)
    hp_diff, hp_i = min_diff(hp_values)

    print(f"x range [{xs[0]:.6g}, {xs[-1]:.6g}], n={len(xs)}")
    print(f"  min h={min(h_values):.12e}")
    print(f"  min h'={min(hp_values):.12e}")
    print(f"  min diff h/x={h_x_diff:.12e} at x={xs[h_x_i]:.6g}")
    print(f"  min diff h'/h={hp_h_diff:.12e} at x={xs[hp_h_i]:.6g}")
    print(f"  min diff x h'/h={x_hp_h_diff:.12e} at x={xs[x_hp_h_i]:.6g}")
    print(f"  min diff h'={hp_diff:.12e} at x={xs[hp_i]:.6g}")


if __name__ == "__main__":
    main()
