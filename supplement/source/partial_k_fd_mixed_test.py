#!/usr/bin/env python3
import argparse

from klm_test import simpson_grid
from partial_k_test import partial_k
from partial_phi_mixed_test import full_partial_h

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_k_fd_mixed_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def fd_h(a, b, omega, nmax, rs, ws, h):
    return (
        partial_k(a + h, b + h, omega, nmax, rs, ws)
        - partial_k(a + h, b - h, omega, nmax, rs, ws)
        - partial_k(a - h, b + h, omega, nmax, rs, ws)
        + partial_k(a - h, b - h, omega, nmax, rs, ws)
    ) / (4 * h * h)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--nmax", type=int, default=2)
    parser.add_argument("--n-grid", type=int, default=41)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--h", type=float, default=0.01)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    for kind in ("analytic", "finite-diff"):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                if kind == "analytic":
                    value = full_partial_h(float(x), float(y), args.omega, args.nmax, rs, ws)
                else:
                    value = fd_h(float(x), float(y), args.omega, args.nmax, rs, ws, args.h)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(
            f"{kind} omega={args.omega:g} nmax={args.nmax} h={args.h:g}: "
            f"min={evals[0]:.12e} max={evals[-1]:.12e}"
        )


if __name__ == "__main__":
    main()
