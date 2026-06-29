#!/usr/bin/env python3
r"""Verify the concrete endpoint Sturm Green identity for the B-model.

Let x=1+tau, c=pi, lambda=c exp(s), and

    F_s(tau) = (exp(-lambda tau)-exp(-c tau))/(lambda-c),
    H_s(tau) = s exp(-lambda tau)/(lambda-c),
    D Y      = (x Y)'.

The endpoint B-kernel is

    K_B(s,t) = int x^(3/2) [
        log(x) D F_s D F_t
        + 1/2 D F_s D H_t
        + 1/2 D H_s D F_t
    ] dtau.

For a coefficient a(tau), define the Sturm expression

    L_a G = -x d/dtau [ a x^(3/2) D G ].

Green's identity gives

    int a x^(3/2) D F D G dtau
      = [a x^(3/2) D G * xF]_{0}^{infty}
        + int F L_a G dtau.

Since every F_s has F_s(0)=0, and log(x)=0 at tau=0, the oriented identity

    K_B(f,u)
      = int F_f L_log F_u dtau
        + 1/2 int F_f L_1 H_u dtau
        + 1/2 int F_u L_1 H_f dtau

has no physical endpoint boundary term.  This is the concrete endpoint Sturm
Green identity needed before the remaining spectral/trace commutation step.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402


def lam_for(c, s):
    return c * mp.e**s


def exp_lam(lam, tau):
    return mp.e ** (-lam * tau)


def f_feature(c, s, tau):
    lam = lam_for(c, s)
    return (exp_lam(lam, tau) - exp_lam(c, tau)) / (lam - c)


def h_feature(c, s, tau):
    lam = lam_for(c, s)
    return s * exp_lam(lam, tau) / (lam - c)


def d_exp(lam, x, tau):
    return (1 - lam * x) * exp_lam(lam, tau)


def d_f_feature(c, s, x, tau):
    lam = lam_for(c, s)
    return (d_exp(lam, x, tau) - d_exp(c, x, tau)) / (lam - c)


def d_h_feature(c, s, x, tau):
    lam = lam_for(c, s)
    return s * d_exp(lam, x, tau) / (lam - c)


def l_exp(lam, x, tau, kind):
    """Return L_a exp(-lambda tau), with a=1 or a=log(x)."""
    if kind == "one":
        a = mp.mpf("1")
        ap = mp.mpf("0")
    elif kind == "log":
        a = mp.log(x)
        ap = 1 / x
    else:
        raise ValueError(kind)

    # M = a x^(3/2) (1-lambda x) exp(-lambda tau)
    # L_a e^{-lambda tau} = -x M'.
    pref = 1 - lam * x
    pref_prime = -lam
    x32 = x ** mp.mpf("1.5")
    x32_prime = mp.mpf("1.5") * x ** mp.mpf("0.5")
    m_no_exp_prime = ap * x32 * pref + a * x32_prime * pref + a * x32 * pref_prime
    m_prime = exp_lam(lam, tau) * (m_no_exp_prime - lam * a * x32 * pref)
    return -x * m_prime


def l_f_feature(c, s, x, tau, kind):
    lam = lam_for(c, s)
    return (l_exp(lam, x, tau, kind) - l_exp(c, x, tau, kind)) / (lam - c)


def l_h_feature(c, s, x, tau, kind):
    lam = lam_for(c, s)
    return s * l_exp(lam, x, tau, kind) / (lam - c)


def combo(c, ss, coeffs, tau, which):
    total = mp.mpf("0")
    for s, coeff in zip(ss, coeffs):
        if which == "F":
            total += coeff * f_feature(c, s, tau)
        elif which == "H":
            total += coeff * h_feature(c, s, tau)
        else:
            raise ValueError(which)
    return total


def combo_d(c, ss, coeffs, x, tau, which):
    total = mp.mpf("0")
    for s, coeff in zip(ss, coeffs):
        if which == "F":
            total += coeff * d_f_feature(c, s, x, tau)
        elif which == "H":
            total += coeff * d_h_feature(c, s, x, tau)
        else:
            raise ValueError(which)
    return total


def combo_l(c, ss, coeffs, x, tau, which, kind):
    total = mp.mpf("0")
    for s, coeff in zip(ss, coeffs):
        if which == "F":
            total += coeff * l_f_feature(c, s, x, tau, kind)
        elif which == "H":
            total += coeff * l_h_feature(c, s, x, tau, kind)
        else:
            raise ValueError(which)
    return total


def direct_density(c, fs, fc, us, uc, r):
    x = mp.e**r
    tau = x - 1
    df = combo_d(c, fs, fc, x, tau, "F")
    du = combo_d(c, us, uc, x, tau, "F")
    dhf = combo_d(c, fs, fc, x, tau, "H")
    dhu = combo_d(c, us, uc, x, tau, "H")
    return x ** mp.mpf("2.5") * (
        r * df * du + mp.mpf("0.5") * (df * dhu + dhf * du)
    )


def green_density(c, fs, fc, us, uc, r):
    x = mp.e**r
    tau = x - 1
    ff = combo(c, fs, fc, tau, "F")
    fu = combo(c, us, uc, tau, "F")
    llog_fu = combo_l(c, us, uc, x, tau, "F", "log")
    l1_hu = combo_l(c, us, uc, x, tau, "H", "one")
    l1_hf = combo_l(c, fs, fc, x, tau, "H", "one")
    # Integrate in r: d tau = x dr.
    return x * (ff * llog_fu + mp.mpf("0.5") * ff * l1_hu + mp.mpf("0.5") * fu * l1_hf)


def boundary_value(c, left_s, left_c, right_s, right_c, r, kind, right_which):
    x = mp.e**r
    tau = x - 1
    if kind == "one":
        a = mp.mpf("1")
    elif kind == "log":
        a = mp.log(x)
    else:
        raise ValueError(kind)
    left_f = combo(c, left_s, left_c, tau, "F")
    right_d = combo_d(c, right_s, right_c, x, tau, right_which)
    return a * x ** mp.mpf("1.5") * right_d * x * left_f


def random_case(args):
    rng = random.Random(args.seed)
    fs = sorted(mp.mpf(args.s_min) + (mp.mpf(args.s_max) - mp.mpf(args.s_min)) * rng.random() for _ in range(args.n))
    us = sorted(mp.mpf(args.s_min) + (mp.mpf(args.s_max) - mp.mpf(args.s_min)) * rng.random() for _ in range(args.n))
    fc = [mp.mpf(rng.uniform(-1, 1)) for _ in range(args.n)]
    uc = [mp.mpf(rng.uniform(-1, 1)) for _ in range(args.n)]
    return fs, fc, us, uc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--s-min", default="0.04")
    parser.add_argument("--s-max", default="1.8")
    parser.add_argument("--rmax", default="14")
    parser.add_argument("--order", type=int, default=50)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    fs, fc, us, uc = random_case(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.rmax), args.order)

    direct = mp.mpf("0")
    green = mp.mpf("0")
    max_density_gap = mp.mpf("0")
    for r, weight in zip(r_nodes, r_weights):
        d = direct_density(c, fs, fc, us, uc, r)
        g = green_density(c, fs, fc, us, uc, r)
        direct += weight * d
        green += weight * g
        max_density_gap = max(max_density_gap, abs(d - g))

    # Boundary terms for the three oriented pieces used in the identity.
    b_log = boundary_value(c, fs, fc, us, uc, 0, "log", "F")
    b_fh = boundary_value(c, fs, fc, us, uc, 0, "one", "H")
    b_uf = boundary_value(c, us, uc, fs, fc, 0, "one", "H")

    print("endpoint Sturm Green identity")
    print(f"  n={args.n} s=[{args.s_min},{args.s_max}] rmax={args.rmax} order={args.order}")
    print("  f nodes:", " ".join(mp.nstr(s, 8) for s in fs))
    print("  u nodes:", " ".join(mp.nstr(s, 8) for s in us))
    print(f"  direct K_B(f,u)       = {mp.nstr(direct, 24)}")
    print(f"  Green oriented total  = {mp.nstr(green, 24)}")
    print(f"  integral defect       = {mp.nstr(direct - green, 12)}")
    print(f"  max density gap       = {mp.nstr(max_density_gap, 12)}")
    print("    (density gap is a total derivative, so it need not be pointwise small)")
    print("  oriented boundary at tau=0:")
    print(f"    log(F_f,F_u)        = {mp.nstr(b_log, 12)}")
    print(f"    one(F_f,H_u)        = {mp.nstr(b_fh, 12)}")
    print(f"    one(F_u,H_f)        = {mp.nstr(b_uf, 12)}")
    print("  identity:")
    print("    K_B(f,u)=int F_f L_log F_u + 1/2 F_f L_1 H_u + 1/2 F_u L_1 H_f")


if __name__ == "__main__":
    main()
