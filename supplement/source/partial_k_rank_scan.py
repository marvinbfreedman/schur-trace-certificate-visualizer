#!/usr/bin/env python3
import argparse

from klm_test import simpson_grid
from partial_k_test import partial_k

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_k_rank_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n-grid", type=int, default=120)
    parser.add_argument("--xmax", type=float, default=4.0)
    parser.add_argument("--nmax", type=int, default=3)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=800)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = partial_k(float(x), float(y), args.omega, args.nmax, rs, ws)
            mat[i, j] = mat[j, i] = value
    evals = np.linalg.eigvalsh(mat)
    max_eval = max(evals[-1], 1.0)
    print(f"omega={args.omega:g} n={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}] K core<={args.nmax}")
    print(f"  min={evals[0]:.12e} max={evals[-1]:.12e}")
    for cutoff in (1e-6, 1e-9, 1e-12, 1e-15):
        print(f"  rank>{cutoff:g}*max(1,lambda_max): {(evals > cutoff * max_eval).sum()}")
    print("  high spectrum:", " ".join(f"{value:.6e}" for value in evals[-20:]))
    print("  low spectrum:", " ".join(f"{value:.6e}" for value in evals[:12]))


if __name__ == "__main__":
    main()
