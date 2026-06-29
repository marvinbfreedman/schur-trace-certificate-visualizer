#!/usr/bin/env python3
import argparse
import math

from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("model_k_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def model_phi(t: float, kind: str) -> float:
    x = t
    if kind == "gaussian":
        return math.exp(-x * x)
    if kind == "hermite2":
        return (1.0 - 2.0 * x * x) * math.exp(-x * x)
    if kind == "hermite4":
        x2 = x * x
        return (1.0 - 4.0 * x2 + 4.0 * x2 * x2 / 3.0) * math.exp(-x2)
    if kind == "two_gauss":
        return math.exp(-x * x) - 0.35 * math.exp(-4.0 * x * x)
    raise ValueError(kind)


def model_k(a: float, b: float, omega: float, kind: str, rs, ws) -> float:
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
            * model_phi(y + u, kind)
            * model_phi(y - u, kind)
        )
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--kinds", default="gaussian,hermite2,hermite4,two_gauss")
    parser.add_argument("--n-grid", type=int, default=81)
    parser.add_argument("--xmax", type=float, default=3.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=800)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for kind in args.kinds.split(","):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = model_k(float(x), float(y), args.omega, kind, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  {kind}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
