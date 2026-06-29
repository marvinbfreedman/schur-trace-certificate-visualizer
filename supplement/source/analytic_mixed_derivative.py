#!/usr/bin/env python3
import argparse
import math

from klm_test import PI, simpson_grid
from mixed_derivative_test import mixed_derivative

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("analytic_mixed_derivative.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def phi_even_derivatives(t: float):
    """Return Phi(t), Phi'(t), Phi''(t) for the even extension used here."""
    sign = 1.0
    x = t
    if x < 0.0:
        sign = -1.0
        x = -x
    if x >= 8.0:
        return 0.0, 0.0, 0.0

    e2 = math.exp(2.0 * x)
    total0 = 0.0
    total1 = 0.0
    total2 = 0.0
    for n in range(1, 80):
        n2 = n * n
        c = PI * n2
        pieces = (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        )
        term_size = 0.0
        for coeff, lam in pieces:
            exp_value = math.exp(lam * x - c * e2)
            log_derivative = lam - 2.0 * c * e2
            second_factor = log_derivative * log_derivative - 4.0 * c * e2
            term0 = coeff * exp_value
            total0 += term0
            total1 += term0 * log_derivative
            total2 += term0 * second_factor
            term_size += abs(term0)
        if n > 4 and term_size < 1e-16 * max(1.0, abs(total0)):
            break

    if t == 0.0:
        return total0, 0.0, total2
    return total0, sign * total1, total2


def w_values(u: float, omega: float):
    ou = omega * u
    cosh = math.cosh(ou)
    sinh = math.sinh(ou)
    w0 = 0.25 * u * cosh
    w1 = 0.25 * (cosh + omega * u * sinh)
    w2 = 0.25 * (2.0 * omega * sinh + omega * omega * u * cosh)
    return w0, w1, w2


def source_st(s: float, t: float, omega: float) -> float:
    fs, fps, _ = phi_even_derivatives(s)
    ft, fpt, _ = phi_even_derivatives(t)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def same_integral(x: float, y: float, omega: float, rs, ws) -> float:
    total = 0.0
    for r, w in zip(rs, ws):
        total += w * source_st(x + r, y + r, omega)
    return total


def reflected_integral(x: float, y: float, omega: float, rs, ws) -> float:
    big = max(x, y)
    small = min(x, y)
    total = 0.0
    for r, w in zip(rs, ws):
        total += w * source_st(big + r, r - small, omega)
    return total


def analytic_h(x: float, y: float, omega: float, rs, ws, parity: str) -> float:
    same = same_integral(x, y, omega, rs, ws)
    reflected = reflected_integral(x, y, omega, rs, ws)
    if parity == "even":
        return same - reflected
    if parity == "odd":
        return same + reflected
    raise ValueError(parity)


def compare(args):
    rs, ws = simpson_grid(args.rmax, args.intervals)
    xs = np.linspace(args.xmin, args.xmax, args.n)
    print(
        f"omega={args.omega:g} n={args.n} x=[{args.xmin:g},{args.xmax:g}] "
        f"fd_h={args.fd_h:g} min_sep={args.min_sep:g}"
    )
    for parity in ("odd", "even"):
        max_abs = 0.0
        max_rel = 0.0
        worst = None
        checked = 0
        for x in xs:
            for y in xs:
                if abs(float(x - y)) < args.min_sep:
                    continue
                analytic = analytic_h(float(x), float(y), args.omega, rs, ws, parity)
                finite = mixed_derivative(float(x), float(y), args.omega, rs, ws, parity, args.fd_h)
                abs_err = abs(analytic - finite)
                rel_err = abs_err / max(1e-16, abs(finite), abs(analytic))
                checked += 1
                if abs_err > max_abs:
                    max_abs = abs_err
                    max_rel = rel_err
                    worst = (float(x), float(y), analytic, finite)
        print(
            f"  {parity}: checked={checked} max_abs={max_abs:.12e} "
            f"max_rel_at_abs={max_rel:.12e}"
        )
        print(
            f"    worst x={worst[0]:.6g} y={worst[1]:.6g} "
            f"analytic={worst[2]:.12e} finite={worst[3]:.12e}"
        )


def spectrum(args):
    rs, ws = simpson_grid(args.rmax, args.intervals)
    xs = np.linspace(args.xmin, args.xmax, args.n)
    print(f"omega={args.omega:g} analytic spectra n={args.n} x=[{args.xmin:g},{args.xmax:g}]")
    for parity in ("odd", "even"):
        mat = np.zeros((args.n, args.n), dtype=float)
        for i, x in enumerate(xs):
            for j, y in enumerate(xs[: i + 1]):
                value = analytic_h(float(x), float(y), args.omega, rs, ws, parity)
                mat[i, j] = mat[j, i] = value
        evals = np.linalg.eigvalsh(mat)
        print(f"  {parity}: min={evals[0]:.12e} max={evals[-1]:.12e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n", type=int, default=20)
    parser.add_argument("--xmin", type=float, default=0.05)
    parser.add_argument("--xmax", type=float, default=2.0)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--fd-h", type=float, default=2e-3)
    parser.add_argument("--min-sep", type=float, default=1e-2)
    parser.add_argument("--mode", choices=("compare", "spectrum"), default="compare")
    args = parser.parse_args()

    if args.mode == "compare":
        compare(args)
    else:
        spectrum(args)


if __name__ == "__main__":
    main()
