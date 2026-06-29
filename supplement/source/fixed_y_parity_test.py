#!/usr/bin/env python3
import argparse

from klm_test import jacobi_eigenvalues_symmetric, phi

try:
    import numpy as np
except ImportError:
    np = None


def min_eig(matrix):
    if np is not None:
        return float(np.linalg.eigvalsh(np.array(matrix, dtype=float))[0])
    return jacobi_eigenvalues_symmetric([row[:] for row in matrix], tol=1e-12)[0]


def fixed_y_parity(y_level: float, n: int, xmax: float):
    xs = [xmax * i / (n - 1) for i in range(n)]
    even = [[0.0 for _ in xs] for _ in xs]
    odd = [[0.0 for _ in xs] for _ in xs]
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            same = 0.0
            if y_level >= 0.5 * (x + y):
                same = phi(y_level + 0.5 * (x - y)) * phi(y_level - 0.5 * (x - y))
            cross = 0.0
            if y_level >= 0.5 * abs(x - y):
                cross = phi(y_level + 0.5 * (x + y)) * phi(y_level - 0.5 * (x + y))
            even[i][j] = even[j][i] = same + cross
            odd[i][j] = odd[j][i] = same - cross
    return min_eig(even), min_eig(odd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--levels", default="0.05,0.1,0.2,0.5,1,2")
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--xmax", type=float, default=8.0)
    args = parser.parse_args()

    for level in [float(x) for x in args.levels.split(",") if x]:
        even, odd = fixed_y_parity(level, args.n, args.xmax)
        print(f"Y={level:g} even_min={even:.12e} odd_min={odd:.12e}")


if __name__ == "__main__":
    main()
