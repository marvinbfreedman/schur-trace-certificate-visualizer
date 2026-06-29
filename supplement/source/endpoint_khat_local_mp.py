#!/usr/bin/env python3
"""High-precision local matrix for Khat=Lhat+Chat.

Also builds the A-branch endpoint kernel

  K_A = int e^(5r/2)[ r AF_i AF_j + 1/2(AF_i AH_j+AH_i AF_j)] dr

and checks the exact transfer identity

  K_B - K_A = 1/2 <F_i,F_j>_{e^(5r/2)dr}.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_bordered_minors_mp import beta, p0  # noqa: E402


def g_value(lam, r):
    x = mp.e**r
    return (1 - lam * x) * mp.e ** (-lam * (x - 1))


def v_value(c, lam, r):
    return (g_value(lam, r) - g_value(c, r)) / (lam - c)


def f_value(c, lam, r):
    x = mp.e**r
    return (mp.e ** (-lam * (x - 1)) - mp.e ** (-c * (x - 1))) / (lam - c)


def q_value(lam, r):
    x = mp.e**r
    return (mp.mpf("1.5") - lam * x) * mp.e ** (-lam * (x - 1))


def af_value(c, lam, r):
    return (q_value(lam, r) - q_value(c, r)) / (lam - c)


def ah_value(c, lam, r):
    return mp.log(lam / c) * q_value(lam, r) / (lam - c)


def l_entry(c, x, y, rmax):
    def integrand(r):
        return (
            mp.e ** (mp.mpf("2.5") * r)
            * r
            * v_value(c, x, r)
            * v_value(c, y, r)
        )

    points = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    z = mp.mpf("0.4")
    while z < rmax:
        points.append(z)
        z *= 2
    points.append(rmax)
    return mp.quad(integrand, points)


def integral_entry(c, x, y, rmax, kind):
    def values(r):
        if kind == "fgram":
            return f_value(c, x, r), f_value(c, y, r), None, None
        if kind == "ka":
            return af_value(c, x, r), af_value(c, y, r), ah_value(c, x, r), ah_value(c, y, r)
        raise ValueError(kind)

    def integrand(r):
        a, b, ax, by = values(r)
        if kind == "fgram":
            return mp.e ** (mp.mpf("2.5") * r) * a * b
        return mp.e ** (mp.mpf("2.5") * r) * (
            r * a * b + mp.mpf("0.5") * (a * by + ax * b)
        )

    points = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    z = mp.mpf("0.4")
    while z < rmax:
        points.append(z)
        z *= 2
    points.append(rmax)
    return mp.quad(integrand, points)


def phat(c, x, y):
    return p0(c, x, y) / ((x - c) * (y - c))


def chat(c, x, y):
    sx = mp.log(x / c)
    sy = mp.log(y / c)
    dx = beta(c, x) / (x - c)
    dy = beta(c, y) / (y - c)
    rx = sx / (x - c)
    ry = sy / (y - c)
    return phat(c, x, y) + mp.mpf("0.5") * (dx * ry + rx * dy)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--rmax", type=str, default="12")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h = mp.mpf(args.h)
    rmax = mp.mpf(args.rmax)
    s = [s0 + i * h for i in range(args.n)]
    lam = [c * mp.e**x for x in s]
    lmat = mp.matrix(args.n)
    cmat = mp.matrix(args.n)
    kmat = mp.matrix(args.n)
    amat = mp.matrix(args.n)
    fgram = mp.matrix(args.n)
    for i in range(args.n):
        for j in range(i + 1):
            le = l_entry(c, lam[i], lam[j], rmax)
            ce = chat(c, lam[i], lam[j])
            ae = integral_entry(c, lam[i], lam[j], rmax, "ka")
            fe = integral_entry(c, lam[i], lam[j], rmax, "fgram")
            lmat[i, j] = lmat[j, i] = le
            cmat[i, j] = cmat[j, i] = ce
            kmat[i, j] = kmat[j, i] = le + ce
            amat[i, j] = amat[j, i] = ae
            fgram[i, j] = fgram[j, i] = fe
    lvals = mp.eigsy(lmat, eigvals_only=True)
    cvals = mp.eigsy(cmat, eigvals_only=True)
    kvals = mp.eigsy(kmat, eigvals_only=True)
    avals = mp.eigsy(amat, eigvals_only=True)
    defect = mp.norm(kmat - amat - mp.mpf("0.5") * fgram)
    print(f"endpoint Khat local mp n={args.n} s0={s0} h={h} dps={args.dps}")
    print("  L eig min/max:", mp.nstr(lvals[0], 20), mp.nstr(lvals[-1], 20))
    print("  A eig min/max:", mp.nstr(avals[0], 20), mp.nstr(avals[-1], 20))
    print("  C eig low:")
    for val in cvals[: min(5, len(cvals))]:
        print(f"    {mp.nstr(val, 25)}")
    print("  K eig low:")
    for val in kvals[: min(8, len(kvals))]:
        print(f"    {mp.nstr(val, 25)}")
    print("  ||K_B-K_A-1/2Fgram||:", mp.nstr(defect, 20))


if __name__ == "__main__":
    main()
