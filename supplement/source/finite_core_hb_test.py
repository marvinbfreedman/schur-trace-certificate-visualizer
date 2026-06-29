#!/usr/bin/env python3
import argparse
import random

import mpmath as mp


def mode_right_derivative(n):
    n2 = n * n
    c = mp.pi * n2
    total = mp.mpf(0)
    for coeff, lam in (
        (4 * mp.pi * mp.pi * n2 * n2, mp.mpf("4.5")),
        (-6 * mp.pi * n2, mp.mpf("2.5")),
    ):
        term = coeff * mp.e ** (-c)
        total += term * (lam - 2 * c)
    return total


def alpha3_zero_slope():
    return -(mode_right_derivative(1) + mode_right_derivative(2)) / mode_right_derivative(3)


def one_sided_transform_term(n, coeff, lam, z):
    c = mp.pi * n * n
    a = (lam + 1j * z) / 2
    return mp.mpf("0.5") * coeff * c ** (-a) * mp.gammainc(a, c, mp.inf)


def xi_finite(z, weights):
    total = 0
    for n, weight in weights.items():
        n2 = n * n
        for coeff, lam in (
            (4 * mp.pi * mp.pi * n2 * n2, mp.mpf("4.5")),
            (-6 * mp.pi * n2, mp.mpf("2.5")),
        ):
            total += weight * (
                one_sided_transform_term(n, coeff, lam, z)
                + one_sided_transform_term(n, coeff, lam, -z)
            )
    return total


def weights_for(kind):
    if kind == "raw3":
        return {1: mp.mpf(1), 2: mp.mpf(1), 3: mp.mpf(1)}
    if kind == "tilde3":
        return {1: mp.mpf(1), 2: mp.mpf(1), 3: alpha3_zero_slope()}
    if kind == "raw2":
        return {1: mp.mpf(1), 2: mp.mpf(1)}
    raise ValueError(kind)


def hb_gap(z, omega, weights):
    plus = abs(xi_finite(z + 1j * omega, weights))
    minus = abs(xi_finite(z - 1j * omega, weights))
    scale = max(plus, minus, mp.mpf("1e-80"))
    return (plus - minus) / scale, plus, minus


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--xmax", default="80")
    parser.add_argument("--ymin", default="0.001")
    parser.add_argument("--ymax", default="4")
    parser.add_argument("--nx", type=int, default=25)
    parser.add_argument("--ny", type=int, default=12)
    parser.add_argument("--random", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    xmax = mp.mpf(args.xmax)
    ymin = mp.mpf(args.ymin)
    ymax = mp.mpf(args.ymax)
    weights = weights_for(args.kind)
    worst = (mp.inf, None, None, None)

    def check(z):
        nonlocal worst
        gap, plus, minus = hb_gap(z, omega, weights)
        if gap < worst[0]:
            worst = (gap, z, plus, minus)

    for i in range(args.nx):
        x = xmax * i / (args.nx - 1)
        for j in range(args.ny):
            y = ymin * (ymax / ymin) ** (mp.mpf(j) / (args.ny - 1))
            check(mp.mpc(x, y))

    random.seed(20260529)
    for _ in range(args.random):
        x = random.random() * float(xmax)
        y = float(ymin) * (float(ymax) / float(ymin)) ** random.random()
        check(mp.mpc(x, y))

    gap, z, plus, minus = worst
    print(f"kind={args.kind} omega={omega} dps={args.dps}")
    print(f"  alpha3={mp.nstr(weights.get(3, mp.mpf(0)), 35)}")
    print(f"  worst_gap={(mp.nstr(gap, 30))}")
    print(f"  at z={mp.nstr(mp.re(z), 18)} + {mp.nstr(mp.im(z), 18)}i")
    print(f"  |Xi(z+iomega)|={mp.nstr(plus, 30)}")
    print(f"  |Xi(z-iomega)|={mp.nstr(minus, 30)}")
    if gap <= 0:
        print("  HB inequality FAILED at sampled point")
    else:
        print("  HB inequality passed sampled points")


if __name__ == "__main__":
    main()
