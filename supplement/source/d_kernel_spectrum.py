#!/usr/bin/env python3
import argparse

from analytic_mixed_derivative import source_st

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("d_kernel_spectrum.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=2.6)
    args = parser.parse_args()

    xs = np.linspace(args.xmin, args.xmax, args.n)
    mat = np.zeros((args.n, args.n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = source_st(float(x), float(y), args.omega)
            mat[i, j] = mat[j, i] = value
    evals = np.linalg.eigvalsh(mat)
    diag = np.diag(mat)
    print(f"omega={args.omega:g} n={args.n} x=[{args.xmin:g},{args.xmax:g}]")
    print(f"  D kernel min={evals[0]:.12e} max={evals[-1]:.12e}")
    print(f"  diagonal min={diag.min():.12e} max={diag.max():.12e}")


if __name__ == "__main__":
    main()
