#!/usr/bin/env python3
import argparse

from partial_phi_mixed_test import partial_source

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_d_kernel_spectrum.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--nmax-values", default="1,2,3,4,8,12")
    parser.add_argument("--n-grid", type=int, default=80)
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=2.6)
    args = parser.parse_args()

    xs = np.linspace(args.xmin, args.xmax, args.n_grid)
    print(f"omega={args.omega:g} D grids n={args.n_grid} x=[{args.xmin:g},{args.xmax:g}]")
    for nmax in (int(text) for text in args.nmax_values.split(",")):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = partial_source(float(x), float(y), args.omega, nmax)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        diag = np.diag(mat)
        print(
            f"  nmax={nmax}: min={evals[0]:.12e} max={evals[-1]:.12e} "
            f"diag_min={diag.min():.12e}"
        )


if __name__ == "__main__":
    main()
