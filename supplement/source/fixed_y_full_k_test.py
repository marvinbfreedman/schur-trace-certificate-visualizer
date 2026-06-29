#!/usr/bin/env python3
import argparse
import math

from partial_k_test import partial_phi

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("fixed_y_full_k_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def fixed_y_kernel(a: float, b: float, y: float, nmax: int) -> float:
    m = 0.5 * (a + b)
    if abs(m) > y:
        return 0.0
    u = 0.5 * (a - b)
    return partial_phi(y + u, nmax) * partial_phi(y - u, nmax)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ys", default="0.2,0.5,1,2,4")
    parser.add_argument("--nmax-values", default="1,2,3")
    parser.add_argument("--n-grid", type=int, default=81)
    parser.add_argument("--xmax", type=float, default=3.0)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    for nmax in (int(text) for text in args.nmax_values.split(",")):
        print(f"nmax={nmax}")
        for y in (float(text) for text in args.ys.split(",")):
            mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
            for i, x in enumerate(xs):
                for j, z in enumerate(xs[: i + 1]):
                    value = fixed_y_kernel(float(x), float(z), y, nmax)
                    mat[i, j] = mat[j, i] = value
            evals = np.linalg.eigvalsh(mat)
            print(f"  y={y:g}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
