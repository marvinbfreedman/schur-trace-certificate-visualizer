#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import w_values
from klm_test import PI, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("scaled_mode3_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def phi_scaled_derivative(t: float, alpha3: float):
    sign = 1.0
    x = t
    if x < 0.0:
        sign = -1.0
        x = -x
    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
    for n in range(1, 4):
        scale = alpha3 if n == 3 else 1.0
        n2 = n * n
        c = PI * n2
        for coeff, lam in (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        ):
            exp_value = math.exp(lam * x - c * e2)
            log_derivative = lam - 2.0 * c * e2
            term0 = scale * coeff * exp_value
            total0 += term0
            total1 += term0 * log_derivative
    if t == 0.0:
        return total0, 0.0
    return total0, sign * total1


def source(s: float, t: float, omega: float, alpha3: float) -> float:
    fs, fps = phi_scaled_derivative(s, alpha3)
    ft, fpt = phi_scaled_derivative(t, alpha3)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def full_h(a: float, b: float, omega: float, alpha3: float, rs, ws) -> float:
    x = abs(a)
    y = abs(b)
    same_sign = a == 0.0 or b == 0.0 or a * b > 0.0
    total = 0.0
    if same_sign:
        for r, w in zip(rs, ws):
            total += w * source(x + r, y + r, omega, alpha3)
    else:
        big = max(x, y)
        small = min(x, y)
        for r, w in zip(rs, ws):
            total += w * source(big + r, r - small, omega, alpha3)
    return total


def right_derivative_at_zero(alpha3: float) -> float:
    _, d = phi_scaled_derivative(1e-300, alpha3)
    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--alphas", default="0,0.25,0.5,0.75,0.9,0.99,1.0,1.01,1.1")
    parser.add_argument("--n-grid", type=int, default=64)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=700)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for alpha in (float(text) for text in args.alphas.split(",")):
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = full_h(float(x), float(y), args.omega, alpha, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(
            f"  alpha3={alpha:.8g}: min={evals[0]:.12e} "
            f"right_deriv0={right_derivative_at_zero(alpha):.12e}"
        )


if __name__ == "__main__":
    main()
