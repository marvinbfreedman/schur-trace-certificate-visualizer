#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import w_values
from klm_test import PI, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("partial_phi_mixed_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def partial_phi_derivative(t: float, nmax: int):
    sign = 1.0
    x = t
    if x < 0.0:
        sign = -1.0
        x = -x
    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
    for n in range(1, nmax + 1):
        n2 = n * n
        c = PI * n2
        for coeff, lam in (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        ):
            exp_value = math.exp(lam * x - c * e2)
            log_derivative = lam - 2.0 * c * e2
            term0 = coeff * exp_value
            total0 += term0
            total1 += term0 * log_derivative
    if t == 0.0:
        return total0, 0.0
    return total0, sign * total1


def partial_source(s: float, t: float, omega: float, nmax: int) -> float:
    fs, fps = partial_phi_derivative(s, nmax)
    ft, fpt = partial_phi_derivative(t, nmax)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def full_partial_h(a: float, b: float, omega: float, nmax: int, rs, ws) -> float:
    x = abs(a)
    y = abs(b)
    same_sign = a == 0.0 or b == 0.0 or a * b > 0.0
    total = 0.0
    if same_sign:
        for r, w in zip(rs, ws):
            total += w * partial_source(x + r, y + r, omega, nmax)
    else:
        big = max(x, y)
        small = min(x, y)
        for r, w in zip(rs, ws):
            total += w * partial_source(big + r, r - small, omega, nmax)
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--nmax-values", default="1,2,3,4,5,8,12")
    parser.add_argument("--n-grid", type=int, default=49)
    parser.add_argument("--xmax", type=float, default=2.4)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=600)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for nmax in (int(text) for text in args.nmax_values.split(",")):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = full_partial_h(float(x), float(y), args.omega, nmax, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  nmax={nmax}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
