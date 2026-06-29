#!/usr/bin/env python3
import argparse

from analytic_mixed_derivative import reflected_integral, same_integral
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("h_component_spectrum.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def fill_matrix(xs, omega: float, rs, ws, kind: str):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            if kind == "same":
                value = same_integral(float(x), float(y), omega, rs, ws)
            elif kind == "reflected":
                value = reflected_integral(float(x), float(y), omega, rs, ws)
            elif kind == "odd":
                value = same_integral(float(x), float(y), omega, rs, ws) + reflected_integral(
                    float(x), float(y), omega, rs, ws
                )
            elif kind == "even":
                value = same_integral(float(x), float(y), omega, rs, ws) - reflected_integral(
                    float(x), float(y), omega, rs, ws
                )
            else:
                raise ValueError(kind)
            mat[i, j] = mat[j, i] = value
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=32)
    parser.add_argument("--xmin", type=float, default=0.03)
    parser.add_argument("--xmax", type=float, default=2.4)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    args = parser.parse_args()

    xs = np.linspace(args.xmin, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} n={args.n} x=[{args.xmin:g},{args.xmax:g}]")
    for kind in ("same", "reflected", "odd", "even"):
        mat = fill_matrix(xs, args.omega, rs, ws, kind)
        evals = np.linalg.eigvalsh(mat)
        print(f"  {kind}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
