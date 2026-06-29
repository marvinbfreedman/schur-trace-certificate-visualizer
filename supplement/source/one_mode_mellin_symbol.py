#!/usr/bin/env python3
import argparse
import cmath
import math


def digamma(z):
    # Shift to a region where the asymptotic expansion is accurate.
    out = 0j
    while abs(z) < 20.0:
        out -= 1.0 / z
        z += 1.0
    inv = 1.0 / z
    inv2 = inv * inv
    # log z - 1/(2z) - sum B_{2k}/(2k z^{2k})
    out += cmath.log(z) - 0.5 * inv
    bernoulli_terms = [
        (1.0 / 6.0, 2),
        (-1.0 / 30.0, 4),
        (1.0 / 42.0, 6),
        (-1.0 / 30.0, 8),
        (5.0 / 66.0, 10),
        (-691.0 / 2730.0, 12),
    ]
    power = inv2
    for bernoulli, order in bernoulli_terms:
        out -= bernoulli * power / order
        power *= inv2
    return out


def log_derivative_symbol(p, tau, c):
    z = complex(p, -tau)
    poly = (2.0 * p - 3.0) ** 2 + 4.0 * tau * tau
    return (
        -2.0 * math.log(c)
        + 2.0 * digamma(z).real
        + 4.0 * (2.0 * p - 3.0) / poly
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=float, default=0.0)
    parser.add_argument("--tau-max", type=float, default=20.0)
    parser.add_argument("--n", type=int, default=2000)
    args = parser.parse_args()

    c = math.pi
    p = args.a + 1.25
    best = (float("inf"), None)
    worst = (-float("inf"), None)
    for i in range(args.n + 1):
        tau = args.tau_max * i / args.n
        value = log_derivative_symbol(p, tau, c)
        if value < best[0]:
            best = (value, tau)
        if value > worst[0]:
            worst = (value, tau)
    print(f"one-mode full Mellin symbol derivative a={args.a:g} p={p:g} c=pi")
    print(f"  min={best[0]:.12e} at tau={best[1]:.6g}")
    print(f"  max={worst[0]:.12e} at tau={worst[1]:.6g}")
    print("  values:")
    for tau in (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0):
        if tau <= args.tau_max:
            print(f"    tau={tau:g}: {log_derivative_symbol(p, tau, c):.12e}")


if __name__ == "__main__":
    main()
