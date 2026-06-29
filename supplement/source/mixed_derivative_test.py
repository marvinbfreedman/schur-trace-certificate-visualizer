#!/usr/bin/env python3
import argparse

from kernel_parity_test import kernel_pm, kernel_pp
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("mixed_derivative_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def sector_kernel(x: float, y: float, omega: float, rs, ws, parity: str) -> float:
    a = kernel_pp(x, y, omega, rs, ws)
    b = kernel_pm(x, y, omega, rs, ws)
    if parity == "odd":
        return a - b
    if parity == "even":
        return a + b
    raise ValueError(parity)


def mixed_derivative(x: float, y: float, omega: float, rs, ws, parity: str, h: float) -> float:
    return (
        sector_kernel(x + h, y + h, omega, rs, ws, parity)
        - sector_kernel(x + h, y - h, omega, rs, ws, parity)
        - sector_kernel(x - h, y + h, omega, rs, ws, parity)
        + sector_kernel(x - h, y - h, omega, rs, ws, parity)
    ) / (4.0 * h * h)


def report(parity: str, xs, omega: float, rs, ws, h: float):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = mixed_derivative(float(x), float(y), omega, rs, ws, parity, h)
            mat[i, j] = mat[j, i] = value
    evals = np.linalg.eigvalsh(mat)
    print(f"  {parity}: mixed min={evals[0]:.12e} max={evals[-1]:.12e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=28)
    parser.add_argument("--xmin", type=float, default=0.03)
    parser.add_argument("--xmax", type=float, default=2.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=500)
    parser.add_argument("--h", type=float, default=2e-3)
    args = parser.parse_args()

    xs = np.linspace(args.xmin, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} n={args.n} x=[{args.xmin:g},{args.xmax:g}] h={args.h:g}")
    report("odd", xs, args.omega, rs, ws, args.h)
    report("even", xs, args.omega, rs, ws, args.h)


if __name__ == "__main__":
    main()
