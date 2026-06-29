#!/usr/bin/env python3
import argparse

from kernel_parity_test import kernel_pm, kernel_pp
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("contraction_spectrum.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmax", type=float, default=8.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--relative-cutoff", type=float, default=1e-12)
    args = parser.parse_args()

    xs = np.linspace(0.0, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    a = np.zeros((args.n, args.n), dtype=float)
    b = np.zeros((args.n, args.n), dtype=float)

    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            av = kernel_pp(float(x), float(y), args.omega, rs, ws)
            bv = kernel_pm(float(x), float(y), args.omega, rs, ws)
            a[i, j] = a[j, i] = av
            b[i, j] = b[j, i] = bv

    evals, basis = np.linalg.eigh(a)
    b_evals = np.linalg.eigvalsh(b)
    even_evals = np.linalg.eigvalsh(a + b)
    odd_evals = np.linalg.eigvalsh(a - b)
    keep = evals > args.relative_cutoff * evals[-1]
    inv_sqrt = basis[:, keep] / np.sqrt(evals[keep])
    contraction = inv_sqrt.T @ b @ inv_sqrt
    spectrum = np.linalg.eigvalsh(contraction)

    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g}")
    print(f"  retained rank={int(keep.sum())} / {args.n}")
    print(f"  A positive min={evals[keep][0]:.12e} max={evals[-1]:.12e}")
    print(f"  B ordinary spectrum min={b_evals[0]:.12e} max={b_evals[-1]:.12e}")
    print(f"  A+B min={even_evals[0]:.12e} A-B min={odd_evals[0]:.12e}")
    print(f"  contraction min={spectrum[0]:.12e} max={spectrum[-1]:.12e}")
    print("  low spectrum:", " ".join(f"{x:.6g}" for x in spectrum[: min(8, len(spectrum))]))
    print("  high spectrum:", " ".join(f"{x:.6g}" for x in spectrum[-min(8, len(spectrum)) :]))


if __name__ == "__main__":
    main()
