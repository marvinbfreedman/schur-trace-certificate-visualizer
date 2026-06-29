#!/usr/bin/env python3
import argparse

from analytic_mixed_derivative import reflected_integral, same_integral
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("h_contraction_spectrum.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=32)
    parser.add_argument("--xmin", type=float, default=0.03)
    parser.add_argument("--xmax", type=float, default=2.4)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--relative-cutoff", type=float, default=1e-12)
    args = parser.parse_args()

    xs = np.linspace(args.xmin, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    c = np.zeros((args.n, args.n), dtype=float)
    rmat = np.zeros((args.n, args.n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            same = same_integral(float(x), float(y), args.omega, rs, ws)
            reflected = reflected_integral(float(x), float(y), args.omega, rs, ws)
            c[i, j] = c[j, i] = same
            rmat[i, j] = rmat[j, i] = reflected

    evals, basis = np.linalg.eigh(c)
    keep = evals > args.relative_cutoff * evals[-1]
    inv_sqrt = basis[:, keep] / np.sqrt(evals[keep])
    contraction = inv_sqrt.T @ rmat @ inv_sqrt
    spectrum = np.linalg.eigvalsh(contraction)
    h_plus = np.linalg.eigvalsh(c - rmat)
    h_minus = np.linalg.eigvalsh(c + rmat)

    print(f"omega={args.omega:g} n={args.n} x=[{args.xmin:g},{args.xmax:g}]")
    print(f"  retained rank={int(keep.sum())} / {args.n}")
    print(f"  C min={evals[keep][0]:.12e} max={evals[-1]:.12e}")
    print(f"  H_+ min={h_plus[0]:.12e} H_- min={h_minus[0]:.12e}")
    print(f"  contraction min={spectrum[0]:.12e} max={spectrum[-1]:.12e}")
    print("  spectrum:", " ".join(f"{value:.6g}" for value in spectrum[: min(12, len(spectrum))]))


if __name__ == "__main__":
    main()
