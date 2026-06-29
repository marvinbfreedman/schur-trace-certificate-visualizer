#!/usr/bin/env python3
import argparse
import mpmath as mp


def w_values(u, omega):
    ou = omega * u
    return (
        mp.mpf("0.25") * u * mp.cosh(ou),
        mp.mpf("0.25") * (mp.cosh(ou) + omega * u * mp.sinh(ou)),
        mp.mpf("0.25") * (2 * omega * mp.sinh(ou) + omega * omega * u * mp.cosh(ou)),
    )


def phi_derivative(t, nmax):
    sign = mp.mpf(1)
    x = t
    if x < 0:
        sign = mp.mpf(-1)
        x = -x
    e2 = mp.e ** (2 * x)
    total0 = mp.mpf(0)
    total1 = mp.mpf(0)
    for n in range(1, nmax + 1):
        n2 = n * n
        c = mp.pi * n2
        for coeff, lam in (
            (4 * mp.pi * mp.pi * n2 * n2, mp.mpf("4.5")),
            (-6 * mp.pi * n2, mp.mpf("2.5")),
        ):
            exp_value = mp.e ** (lam * x - c * e2)
            log_derivative = lam - 2 * c * e2
            term0 = coeff * exp_value
            total0 += term0
            total1 += term0 * log_derivative
    if t == 0:
        return total0, mp.mpf(0)
    return total0, sign * total1


def source(s, t, omega, nmax):
    fs, fps = phi_derivative(s, nmax)
    ft, fpt = phi_derivative(t, nmax)
    w0, w1, w2 = w_values(s + t, omega)
    return w2 * fs * ft + w1 * (fps * ft + fs * fpt) + w0 * fps * fpt


def simpson_grid(rmax, intervals):
    if intervals % 2:
        raise ValueError("intervals must be even")
    h = rmax / intervals
    rs = [h * i for i in range(intervals + 1)]
    ws = []
    for i in range(intervals + 1):
        if i == 0 or i == intervals:
            factor = 1
        elif i % 2:
            factor = 4
        else:
            factor = 2
        ws.append(h * factor / 3)
    return rs, ws


def full_h(a, b, omega, nmax, rs, ws):
    x = abs(a)
    y = abs(b)
    same_sign = a == 0 or b == 0 or a * b > 0
    total = mp.mpf(0)
    if same_sign:
        for r, w in zip(rs, ws):
            total += w * source(x + r, y + r, omega, nmax)
    else:
        big = max(x, y)
        small = min(x, y)
        for r, w in zip(rs, ws):
            total += w * source(big + r, r - small, omega, nmax)
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--n-grid", type=int, default=18)
    parser.add_argument("--xmax", default="2.6")
    parser.add_argument("--nmax", type=int, default=3)
    parser.add_argument("--rmax", default="9.0")
    parser.add_argument("--intervals", type=int, default=320)
    parser.add_argument("--dps", type=int, default=70)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    xmax = mp.mpf(args.xmax)
    rmax = mp.mpf(args.rmax)
    rs, ws = simpson_grid(rmax, args.intervals)
    xs = [-xmax + 2 * xmax * i / (args.n_grid - 1) for i in range(args.n_grid)]
    mat = mp.matrix(args.n_grid)
    for i, x in enumerate(xs):
        for j in range(i + 1):
            value = full_h(xs[i], xs[j], omega, args.nmax, rs, ws)
            mat[i, j] = value
            mat[j, i] = value
    evals = mp.eigsy(mat, eigvals_only=True)
    print(f"omega={omega} n={args.n_grid} x=[{-xmax},{xmax}] nmax={args.nmax} dps={args.dps}")
    print("  min=", mp.nstr(evals[0], 30))
    print("  max=", mp.nstr(evals[-1], 30))
    print("  low=", " ".join(mp.nstr(value, 12) for value in evals[: min(8, len(evals))]))


if __name__ == "__main__":
    main()
