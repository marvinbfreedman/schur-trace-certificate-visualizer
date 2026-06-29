#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import phi_even_derivatives
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("h_lg_component_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def g_derivatives(t: float, omega: float, sign: int):
    f, fp, _ = phi_even_derivatives(t)
    exp_value = math.exp(sign * omega * t)
    g = exp_value * f
    gp = exp_value * (sign * omega * f + fp)
    return g, gp


def c_lg_integrand(s: float, t: float, omega: float, sign: int) -> float:
    gs, gps = g_derivatives(s, omega, sign)
    gt, gpt = g_derivatives(t, omega, sign)
    return 0.5 * gps * gt + 0.5 * gs * gpt + 0.5 * (s + t) * gps * gpt


def kernel_component(x: float, y: float, omega: float, sign: int, rs, ws) -> float:
    total = 0.0
    for r, w in zip(rs, ws):
        total += w * c_lg_integrand(x + r, y + r, omega, sign)
    return total


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
    mats = []
    for sign, label in ((1, "g_+"), (-1, "g_-")):
        mat = np.zeros((args.n, args.n), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = kernel_component(float(x), float(y), args.omega, sign, rs, ws)
                mat[i, j] = mat[j, i] = value
        mats.append(mat)
        evals = np.linalg.eigvalsh(mat)
        print(f"  {label}: min={evals[0]:.12e} max={evals[-1]:.12e}")
    combined = 0.25 * (mats[0] + mats[1])
    evals = np.linalg.eigvalsh(combined)
    print(f"  C from components: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
