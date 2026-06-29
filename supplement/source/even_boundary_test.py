#!/usr/bin/env python3
import argparse

from kernel_parity_test import kernel_pm, kernel_pp
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("even_boundary_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def even_kernel(x: float, y: float, omega: float, rs, ws) -> float:
    return kernel_pp(x, y, omega, rs, ws) + kernel_pm(x, y, omega, rs, ws)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--xmax", type=float, default=2.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    args = parser.parse_args()

    xs = np.linspace(0.0, args.xmax, args.n)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    f = np.array([even_kernel(float(x), 0.0, args.omega, rs, ws) for x in xs])
    c = even_kernel(0.0, 0.0, args.omega, rs, ws)
    boundary = f[:, None] + f[None, :] - c
    evals = np.linalg.eigvalsh(boundary)
    print(f"omega={args.omega:g} n={args.n} xmax={args.xmax:g}")
    print(f"  E(0,0)={c:.12e}")
    print(f"  boundary min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
