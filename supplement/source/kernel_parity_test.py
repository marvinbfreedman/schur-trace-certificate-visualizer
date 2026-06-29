#!/usr/bin/env python3
import argparse
import math
import time

from klm_test import jacobi_eigenvalues_symmetric, phi, simpson_grid

try:
    import numpy as np
except ImportError:
    np = None


def kernel_pp(x: float, y: float, omega: float, rs, ws) -> float:
    half_sum = 0.5 * (x + y)
    total = 0.0
    for r, w in zip(rs, ws):
        z = r + half_sum
        total += w * z * math.cosh(2.0 * omega * z) * phi(r + x) * phi(r + y)
    return 0.5 * total


def kernel_pm(x: float, y: float, omega: float, rs, ws) -> float:
    big = max(x, y)
    small = min(x, y)
    half_diff = 0.5 * (big - small)
    total = 0.0
    for r, w in zip(rs, ws):
        z = r + half_diff
        total += w * z * math.cosh(2.0 * omega * z) * phi(r + big) * phi(r - small)
    return 0.5 * total


def min_eig_real_symmetric(matrix):
    if np is not None:
        return float(np.linalg.eigvalsh(np.array(matrix, dtype=float))[0])
    return jacobi_eigenvalues_symmetric([row[:] for row in matrix], tol=1e-12)[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmax", type=float, default=8.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=900)
    args = parser.parse_args()

    start = time.time()
    xs = [args.xmax * i / (args.n - 1) for i in range(args.n)]
    rs, ws = simpson_grid(args.rmax, args.intervals)
    even = [[0.0 for _ in xs] for _ in xs]
    odd = [[0.0 for _ in xs] for _ in xs]

    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            pp = kernel_pp(x, y, args.omega, rs, ws)
            pm = kernel_pm(x, y, args.omega, rs, ws)
            even_value = pp + pm
            odd_value = pp - pm
            even[i][j] = even_value
            even[j][i] = even_value
            odd[i][j] = odd_value
            odd[j][i] = odd_value

    backend = "numpy" if np is not None else "pure-python"
    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g} backend={backend}")
    print(f"  even lambda_min={min_eig_real_symmetric(even):.12e}")
    print(f"  odd  lambda_min={min_eig_real_symmetric(odd):.12e}")
    print(f"  elapsed={time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
