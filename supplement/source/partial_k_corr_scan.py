#!/usr/bin/env python3
import argparse
import random

from klm_test import simpson_grid
from partial_k_test import partial_k

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_k_corr_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def corr_spectrum(xs, omega, nmax, rs, ws):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    diag = np.zeros(n, dtype=float)
    for i, x in enumerate(xs):
        diag[i] = partial_k(float(x), float(x), omega, nmax, rs, ws)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            denom = (diag[i] * diag[j]) ** 0.5
            if denom == 0.0:
                value = 0.0
            else:
                value = partial_k(float(x), float(y), omega, nmax, rs, ws) / denom
            mat[i, j] = mat[j, i] = value
    return np.linalg.eigvalsh(mat), diag


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--nmax", type=int, default=3)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--points", type=int, default=12)
    parser.add_argument("--xmax", type=float, default=20.0)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=900)
    args = parser.parse_args()

    rs, ws = simpson_grid(args.rmax, args.intervals)
    random.seed(20260529)
    worst = (float("inf"), None, None)
    for trial in range(args.trials):
        if trial == 0:
            xs = np.linspace(-args.xmax, args.xmax, args.points)
        elif trial == 1:
            xs = np.r_[
                np.linspace(-args.xmax, -0.01, args.points // 2),
                np.linspace(0.01, args.xmax, args.points - args.points // 2),
            ]
        else:
            xs = np.array(sorted(random.uniform(-args.xmax, args.xmax) for _ in range(args.points)))
        evals, diag = corr_spectrum(xs, args.omega, args.nmax, rs, ws)
        if evals[0] < worst[0]:
            worst = (float(evals[0]), xs.copy(), diag.copy())
    print(
        f"omega={args.omega:g} nmax={args.nmax} trials={args.trials} "
        f"points={args.points} x=[{-args.xmax:g},{args.xmax:g}]"
    )
    print(f"  worst corr min={worst[0]:.12e}")
    print("  xs:", " ".join(f"{x:.6g}" for x in worst[1]))
    print("  diag range:", f"{worst[2].min():.12e}", f"{worst[2].max():.12e}")


if __name__ == "__main__":
    main()
