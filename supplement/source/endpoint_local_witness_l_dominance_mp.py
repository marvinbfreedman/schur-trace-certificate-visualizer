#!/usr/bin/env python3
"""Evaluate L domination on the local high-precision P0 row-null witness."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_bordered_minors_mp import beta, p0  # noqa: E402
from endpoint_row_null_local_mp import null_basis, rows  # noqa: E402


def g_value(lam, r):
    x = mp.e**r
    return (1 - lam * x) * mp.e ** (-lam * (x - 1))


def v_value(c, lam, r):
    return (g_value(lam, r) - g_value(c, r)) / (lam - c)


def khat(c, x, y):
    return p0(c, x, y) / ((x - c) * (y - c))


def l_quadratic(c, lam, alpha, rmax):
    def integrand(r):
        v = mp.fsum(alpha[i] * v_value(c, lam[i], r) for i in range(len(lam)))
        return mp.e ** (mp.mpf("2.5") * r) * r * v * v

    points = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    x = mp.mpf("0.4")
    while x < rmax:
        points.append(x)
        x *= 2
    points.append(rmax)
    return mp.quad(integrand, points)


def quadratic(c, lam, alpha):
    total = mp.mpf("0")
    for i in range(len(lam)):
        for j in range(len(lam)):
            total += alpha[i] * alpha[j] * khat(c, lam[i], lam[j])
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=120)
    parser.add_argument("--rmax", type=str, default="12")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h = mp.mpf(args.h)
    s = [s0 + i * h for i in range(args.n)]
    lam = [c * mp.e**x for x in s]
    e, d = rows(c, lam)
    basis = null_basis(e, d)
    m = len(basis)
    q = mp.matrix(m)
    for a in range(m):
        for b in range(m):
            total = mp.mpf("0")
            for i in range(args.n):
                for j in range(args.n):
                    total += basis[a][i] * basis[b][j] * khat(c, lam[i], lam[j])
            q[a, b] = total
    eigvals, eigvecs = mp.eigsy(q, eigvals_only=False)
    idx = min(range(len(eigvals)), key=lambda k: eigvals[k])
    alpha = [mp.mpf("0")] * args.n
    for a in range(m):
        for i in range(args.n):
            alpha[i] += basis[a][i] * eigvecs[a, idx]

    pval = quadratic(c, lam, alpha)
    lval = l_quadratic(c, lam, alpha, mp.mpf(args.rmax))
    print(
        f"endpoint local witness L dominance n={args.n} s0={s0} h={h} dps={args.dps}"
    )
    print(f"  P0hat witness = {mp.nstr(pval, 25)}")
    print(f"  Lhat witness  = {mp.nstr(lval, 25)}")
    print(f"  L/(-P0)       = {mp.nstr(lval/(-pval), 25)}")
    print(f"  K=L+P0        = {mp.nstr(lval+pval, 25)}")


if __name__ == "__main__":
    main()
