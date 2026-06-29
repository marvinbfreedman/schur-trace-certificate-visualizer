#!/usr/bin/env python3
import argparse

from klm_test import simpson_grid
from wedge_scalar_test import anti_diagonal_derivative, reflected_source

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("slice_kernel_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def odd_density(c: float, v: float, omega: float, rs, ws, step: float) -> float:
    return anti_diagonal_derivative(c, v, omega, rs, ws, step) - reflected_source(c, v, omega)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--v-values", default="0,0.05,0.1,0.2,0.4,0.8")
    parser.add_argument("--n", type=int, default=30)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=500)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g}")
    for text in args.v_values.split(","):
        v = float(text)
        if v >= args.xmax:
            continue
        xs = np.linspace(v, args.xmax, args.n)
        mat = np.zeros((args.n, args.n), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = odd_density(float(x + y), v, args.omega, rs, ws, args.h)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  v={v:.6g}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
