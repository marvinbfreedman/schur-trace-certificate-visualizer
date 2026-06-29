#!/usr/bin/env python3
import argparse
import math

from analytic_mixed_derivative import w_values
from klm_test import PI, simpson_grid
from second_order_cosh_ibp import finite_modes, mode_right_derivative

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_defect_decomp.py requires numpy; run with python, not python3") from exc


def positive_values(t, modes):
    e2 = math.exp(2.0 * t)
    total = 0.0
    total_p = 0.0
    for n, weight in modes:
        n2 = n * n
        c = PI * n2
        for coeff, lam in (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        ):
            term = weight * coeff * math.exp(lam * t - c * e2)
            total += term
            total_p += term * (lam - 2.0 * c * e2)
    return total, total_p


def even_values(t, modes):
    if t > 0.0:
        return positive_values(t, modes)
    if t < 0.0:
        value, deriv = positive_values(-t, modes)
        return value, -deriv
    value, _ = positive_values(0.0, modes)
    return value, 0.0


def right_slope(modes):
    return sum(weight * mode_right_derivative(n) for n, weight in modes)


def source(s, t, omega, modes):
    fs, fps = even_values(s, modes)
    ft, fpt = even_values(t, modes)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def sign0(value):
    if value > 0.0:
        return 1.0
    if value < 0.0:
        return -1.0
    return 0.0


def defect_amplitude(s, t, omega, modes):
    # This is the coefficient of the jump in the crossing variable t:
    # W'(s+t) F(s) + W(s+t) F'(s), with s > 0 in the reflected kernel.
    fs, fps = positive_values(s, modes)
    w0, w1, _ = w_values(s + t, omega)
    return w1 * fs + w0 * fps


def residual_source(s, t, omega, modes, slope):
    return source(s, t, omega, modes) - slope * sign0(t) * defect_amplitude(
        s, t, omega, modes
    )


def reflected_decomp(x, y, omega, modes, rs, ws):
    big = max(abs(x), abs(y))
    small = min(abs(x), abs(y))
    slope = right_slope(modes)
    direct = 0.0
    residual = 0.0
    defect = 0.0
    for r, w in zip(rs, ws):
        s = big + r
        t = r - small
        direct += w * source(s, t, omega, modes)
        residual += w * residual_source(s, t, omega, modes, slope)
        defect += w * sign0(t) * defect_amplitude(s, t, omega, modes)
    return direct, residual, defect


def full_decomp(a, b, omega, modes, rs, ws):
    same_sign = a == 0.0 or b == 0.0 or a * b > 0.0
    slope = right_slope(modes)
    if same_sign:
        x = abs(a)
        y = abs(b)
        direct = 0.0
        for r, w in zip(rs, ws):
            direct += w * source(x + r, y + r, omega, modes)
        return direct, direct, 0.0
    direct, residual, defect = reflected_decomp(a, b, omega, modes, rs, ws)
    return direct, residual, defect


def jump_check(x, y, omega, modes, eps_values):
    big = max(abs(x), abs(y))
    small = min(abs(x), abs(y))
    s0 = big + small
    slope = right_slope(modes)
    amp0 = defect_amplitude(s0, 0.0, omega, modes)
    predicted = 2.0 * slope * amp0
    rows = []
    for eps in eps_values:
        plus = source(s0, eps, omega, modes)
        minus = source(s0, -eps, omega, modes)
        res_plus = residual_source(s0, eps, omega, modes, slope)
        res_minus = residual_source(s0, -eps, omega, modes, slope)
        rows.append(
            {
                "eps": eps,
                "source_jump": plus - minus,
                "predicted": predicted,
                "residual_jump": res_plus - res_minus,
            }
        )
    return rows


def matrix_report(args, modes, rs, ws):
    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    direct = np.zeros((args.n_grid, args.n_grid), dtype=float)
    residual = np.zeros_like(direct)
    defect = np.zeros_like(direct)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            dval, rval, jval = full_decomp(float(x), float(y), args.omega, modes, rs, ws)
            direct[i, j] = direct[j, i] = dval
            residual[i, j] = residual[j, i] = rval
            defect[i, j] = defect[j, i] = jval

    slope = right_slope(modes)
    reconstructed = residual + slope * defect
    err = direct - reconstructed
    direct_vals, direct_vecs = np.linalg.eigh((direct + direct.T) / 2.0)
    residual_vals = np.linalg.eigvalsh((residual + residual.T) / 2.0)
    defect_vals = np.linalg.eigvalsh((slope * defect + (slope * defect).T) / 2.0)
    witness = direct_vecs[:, 0]
    print(f"matrix grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    print(
        f"  identity max_abs={np.max(np.abs(err)):.12e} "
        f"rel={np.max(np.abs(err)) / max(1e-300, np.max(np.abs(direct))):.12e}"
    )
    print(
        f"  eig direct min={direct_vals[0]:.12e} "
        f"residual min={residual_vals[0]:.12e} "
        f"slope*defect min={defect_vals[0]:.12e}"
    )
    print(
        f"  direct-min witness qforms: direct={float(witness @ direct @ witness):.12e} "
        f"residual={float(witness @ residual @ witness):.12e} "
        f"slope*defect={float(witness @ (slope * defect) @ witness):.12e}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="raw2")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--x", type=float, default=0.7)
    parser.add_argument("--y", type=float, default=0.2)
    parser.add_argument("--eps", default="1e-2,1e-3,1e-4,1e-5")
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=900)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n-grid", type=int, default=48)
    args = parser.parse_args()

    modes = finite_modes(args.kind)
    slope = right_slope(modes)
    print(
        f"endpoint defect kind={args.kind} omega={args.omega:g} "
        f"slope F'(0+)={slope:.12e}"
    )
    for row in jump_check(
        args.x, args.y, args.omega, modes, [float(text) for text in args.eps.split(",")]
    ):
        print(
            f"  eps={row['eps']:.1e} source_jump={row['source_jump']:.12e} "
            f"pred={row['predicted']:.12e} "
            f"residual_jump={row['residual_jump']:.12e}"
        )

    rs, ws = simpson_grid(args.rmax, args.intervals)
    direct, residual, defect = reflected_decomp(args.x, args.y, args.omega, modes, rs, ws)
    print(
        f"reflected point x={args.x:g} y={args.y:g}: "
        f"direct={direct:.12e} residual={residual:.12e} "
        f"slope*defect={slope * defect:.12e} "
        f"err={direct - residual - slope * defect:.12e}"
    )
    matrix_report(args, modes, rs, ws)


if __name__ == "__main__":
    main()
