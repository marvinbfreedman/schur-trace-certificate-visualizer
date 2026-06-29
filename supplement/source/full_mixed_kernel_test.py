#!/usr/bin/env python3
import argparse

from analytic_mixed_derivative import reflected_integral, same_integral
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("full_mixed_kernel_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def full_h(a: float, b: float, omega: float, rs, ws) -> float:
    if a == 0.0 or b == 0.0 or a * b > 0.0:
        return same_integral(abs(a), abs(b), omega, rs, ws)
    return reflected_integral(abs(a), abs(b), omega, rs, ws)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=65)
    parser.add_argument("--xmax", type=float, default=2.4)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    mat = np.zeros((args.n, args.n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = full_h(float(x), float(y), args.omega, rs, ws)
            mat[i, j] = mat[j, i] = value
    evals = np.linalg.eigvalsh(mat)
    print(f"omega={args.omega:g} n={args.n} x=[{-args.xmax:g},{args.xmax:g}]")
    print(f"  full mixed kernel min={evals[0]:.12e} max={evals[-1]:.12e}")
    print("  low spectrum:", " ".join(f"{value:.6g}" for value in evals[: min(10, len(evals))]))


if __name__ == "__main__":
    main()
