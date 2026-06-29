#!/usr/bin/env python3
import argparse

from core_witness_test import matrix_for
from klm_test import simpson_grid
from term_pair_mixed_test import full_sym_pair_h

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("pair_witness_decomp.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def pair_matrix(xs, omega, n, m, rs, ws):
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = full_sym_pair_h(float(x), float(y), omega, n, m, rs, ws)
            mat[i, j] = mat[j, i] = value
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n-grid", type=int, default=64)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=1000)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    core2 = matrix_for(xs, args.omega, rs, ws, "core2")
    evals, basis = np.linalg.eigh(core2)
    witness = basis[:, 0]

    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    print(f"  lambda_min core2={evals[0]:.12e}")
    for n, m in ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3)):
        mat = pair_matrix(xs, args.omega, n, m, rs, ws)
        q = float(witness @ mat @ witness)
        print(f"  witness pair ({n},{m})={q:.12e}")


if __name__ == "__main__":
    main()
