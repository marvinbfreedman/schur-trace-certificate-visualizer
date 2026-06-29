#!/usr/bin/env python3
import argparse
import math

from klm_test import PI, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_k_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def partial_phi(t: float, nmax: int) -> float:
    x = abs(t)
    e2 = math.exp(2.0 * x)
    total = 0.0
    for n in range(1, nmax + 1):
        n2 = n * n
        c = PI * n2
        total += (
            4.0 * PI * PI * n2 * n2 * math.exp(4.5 * x - c * e2)
            - 6.0 * PI * n2 * math.exp(2.5 * x - c * e2)
        )
    return total


def partial_k(a: float, b: float, omega: float, nmax: int, rs, ws) -> float:
    m = 0.5 * (a + b)
    u = 0.5 * (a - b)
    lower = abs(m)
    total = 0.0
    for r, w in zip(rs, ws):
        y = lower + r
        total += (
            0.5
            * w
            * y
            * math.cosh(2.0 * omega * y)
            * partial_phi(y + u, nmax)
            * partial_phi(y - u, nmax)
        )
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--nmax-values", default="1,2,3,4,8")
    parser.add_argument("--n-grid", type=int, default=49)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=700)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for nmax in (int(text) for text in args.nmax_values.split(",")):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = partial_k(float(x), float(y), args.omega, nmax, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  nmax={nmax}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
