#!/usr/bin/env python3
import argparse
import math

from klm_test import phi, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("lg_component_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def g_value(t: float, omega: float, sign: int) -> float:
    return math.exp(sign * omega * t) * phi(t)


def kernel_lg(x: float, y: float, omega: float, sign: int, rs, ws) -> float:
    total = 0.0
    for r, w in zip(rs, ws):
        sx = r + x
        sy = r + y
        total += w * 0.5 * (sx + sy) * g_value(sx, omega, sign) * g_value(sy, omega, sign)
    return total


def matrix_lg(xs, omega: float, sign: int, rs, ws):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = kernel_lg(float(x), float(y), omega, sign, rs, ws)
            mat[i, j] = mat[j, i] = value
    return mat


def report_component(xs, omega: float, sign: int, rs, ws):
    mat = matrix_lg(xs, omega, sign, rs, ws)
    evals = np.linalg.eigvalsh(mat)
    label = "g_+" if sign > 0 else "g_-"
    print(
        f"  {label}: lambda_min={evals[0]:.12e} "
        f"lambda_max={evals[-1]:.12e} cond_floor={evals[0] / evals[-1]:.12e}"
    )
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmax", type=float, default=8.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=900)
    args = parser.parse_args()

    xs = np.linspace(0.0, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)

    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g}")
    plus = report_component(xs, args.omega, 1, rs, ws)
    minus = report_component(xs, args.omega, -1, rs, ws)
    combined = 0.25 * (plus + minus)
    evals = np.linalg.eigvalsh(combined)
    print(f"  A from components: lambda_min={evals[0]:.12e} lambda_max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
