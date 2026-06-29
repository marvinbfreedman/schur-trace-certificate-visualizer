#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import w_values
from klm_test import PI, simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("ranged_phi_mixed_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def ranged_phi_derivative(t: float, nmin: int, nmax: int):
    sign = 1.0
    x = t
    if x < 0.0:
        sign = -1.0
        x = -x
    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
    for n in range(nmin, nmax + 1):
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


def ranged_source(s: float, t: float, omega: float, nmin: int, nmax: int) -> float:
    fs, fps = ranged_phi_derivative(s, nmin, nmax)
    ft, fpt = ranged_phi_derivative(t, nmin, nmax)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def full_ranged_h(a: float, b: float, omega: float, nmin: int, nmax: int, rs, ws) -> float:
    x = abs(a)
    y = abs(b)
    same_sign = a == 0.0 or b == 0.0 or a * b > 0.0
    total = 0.0
    if same_sign:
        for r, w in zip(rs, ws):
            total += w * ranged_source(x + r, y + r, omega, nmin, nmax)
    else:
        big = max(x, y)
        small = min(x, y)
        for r, w in zip(rs, ws):
            total += w * ranged_source(big + r, r - small, omega, nmin, nmax)
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--ranges", default="1:1,2:2,1:2,2:8,3:8,1:8")
    parser.add_argument("--n-grid", type=int, default=49)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=600)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    for spec in args.ranges.split(","):
        nmin_text, nmax_text = spec.split(":")
        nmin = int(nmin_text)
        nmax = int(nmax_text)
        mat = np.zeros((args.n_grid, args.n_grid), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = full_ranged_h(float(x), float(y), args.omega, nmin, nmax, rs, ws)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  range {nmin}:{nmax}: min={evals[0]:.12e} max={evals[-1]:.12e}")


if __name__ == "__main__":
    main()
