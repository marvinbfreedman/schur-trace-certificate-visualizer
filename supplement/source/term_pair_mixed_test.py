#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import w_values
from klm_test import PI, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("term_pair_mixed_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def phi_term_derivative(n: int, t: float):
    sign = 1.0
    x = t
    if x < 0.0:
        sign = -1.0
        x = -x
    n2 = n * n
    c = PI * n2
    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
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


def source_pair(s: float, t: float, omega: float, n: int, m: int) -> float:
    fs, fps = phi_term_derivative(n, s)
    gt, gpt = phi_term_derivative(m, t)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * gt + w1 * (fps * gt + fs * gpt) + w0 * fps * gpt


def full_pair_h(a: float, b: float, omega: float, n: int, m: int, rs, ws) -> float:
    same_sign = a == 0.0 or b == 0.0 or a * b > 0.0
    x = abs(a)
    y = abs(b)
    if same_sign:
        total = 0.0
        for r, w in zip(rs, ws):
            total += w * source_pair(x + r, y + r, omega, n, m)
        return total

    big = max(x, y)
    small = min(x, y)
    total = 0.0
    for r, w in zip(rs, ws):
        total += w * source_pair(big + r, r - small, omega, n, m)
    return total


def full_sym_pair_h(a: float, b: float, omega: float, n: int, m: int, rs, ws) -> float:
    if n == m:
        return full_pair_h(a, b, omega, n, m, rs, ws)
    return full_pair_h(a, b, omega, n, m, rs, ws) + full_pair_h(a, b, omega, m, n, rs, ws)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--pairs", default="1,1;1,2;2,2;1,3;2,3;3,3")
    parser.add_argument("--n-grid", type=int, default=49)
    parser.add_argument("--xmax", type=float, default=2.4)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=600)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for pair in args.pairs.split(";"):
        n_text, m_text = pair.split(",")
        n = int(n_text)
        m = int(m_text)
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = full_sym_pair_h(float(x), float(y), args.omega, n, m, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  pair ({n},{m}): min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
