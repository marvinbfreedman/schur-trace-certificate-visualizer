#!/usr/bin/env python3
"""Confluent Taylor-jet tests for the endpoint A/B transfer kernel.

For an analytic positive kernel K(s,t), every confluent jet matrix

  [ d_s^i d_t^j K(s0,t0) / (i! j!) ]_{i,j=0}^{n-1},  t0=s0,

must be positive semidefinite.  This script builds that matrix at the
integrand level by formal Taylor series in delta=s-s0.  This avoids both
close-node conditioning and slow high-order numerical differentiation.

The endpoint integrands are, with lambda=c exp(s), x=exp(r), tau=x-1,

  F = (exp(-lambda tau)-exp(-c tau))/(lambda-c),
  V_B = ((1-lambda x)exp(-lambda tau)-(1-cx)exp(-c tau))/(lambda-c),
  W_B = s (1-lambda x)exp(-lambda tau)/(lambda-c),

and similarly

  V_A = ((3/2-lambda x)exp(-lambda tau)
         -(3/2-cx)exp(-c tau))/(lambda-c),
  W_A = s (3/2-lambda x)exp(-lambda tau)/(lambda-c).

The matrices are

  FGram = int exp(5r/2) F_i F_j dr,
  K_A   = int exp(5r/2) [r V_Ai V_Aj + 1/2(V_Ai W_Aj+W_Ai V_Aj)] dr,
  K_B   = int exp(5r/2) [r V_Bi V_Bj + 1/2(V_Bi W_Bj+W_Bi V_Bj)] dr.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


def zero(n):
    return [mp.mpf("0") for _ in range(n)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [k * x for x in a]


def mul(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = mp.fsum(a[j] * b[i - j] for j in range(i + 1))
    return out


def div(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = (
            a[i] - mp.fsum(b[j] * out[i - j] for j in range(1, i + 1))
        ) / b[0]
    return out


def exp_series(a):
    n = len(a)
    out = zero(n)
    out[0] = mp.e**a[0]
    for i in range(1, n):
        out[i] = mp.fsum(k * a[k] * out[i - k] for k in range(1, i + 1)) / i
    return out


def lambda_series(c, s0, n):
    lam0 = c * mp.e**s0
    return [lam0 / mp.factorial(i) for i in range(n)]


def endpoint_series(c, s0, r, n):
    x = mp.e**r
    tau = x - 1
    lam = lambda_series(c, s0, n)
    denom = lam[:]
    denom[0] -= c

    exp_lam = exp_series(scale(lam, -tau))
    exp_c = mp.e ** (-c * tau)

    one = zero(n)
    one[0] = 1
    sser = zero(n)
    sser[0] = s0
    if n > 1:
        sser[1] = 1

    f_num = exp_lam[:]
    f_num[0] -= exp_c
    f = div(f_num, denom)

    gb_pref = sub(one, scale(lam, x))
    gb = mul(gb_pref, exp_lam)
    gb_c = (1 - c * x) * exp_c
    vb_num = gb[:]
    vb_num[0] -= gb_c
    vb = div(vb_num, denom)
    wb = div(mul(sser, gb), denom)

    ga_pref = sub(scale(one, mp.mpf("1.5")), scale(lam, x))
    ga = mul(ga_pref, exp_lam)
    ga_c = (mp.mpf("1.5") - c * x) * exp_c
    va_num = ga[:]
    va_num[0] -= ga_c
    va = div(va_num, denom)
    wa = div(mul(sser, ga), denom)

    return f, va, wa, vb, wb


def contribution(kind, c, s0, r, n):
    f, va, wa, vb, wb = endpoint_series(c, s0, r, n)
    mat = mp.matrix(n)
    weight = mp.e ** (mp.mpf("2.5") * r)
    for i in range(n):
        for j in range(i + 1):
            if kind == "fgram":
                val = f[i] * f[j]
            elif kind == "ka":
                val = r * va[i] * va[j] + mp.mpf("0.5") * (
                    va[i] * wa[j] + wa[i] * va[j]
                )
            elif kind == "kb":
                val = r * vb[i] * vb[j] + mp.mpf("0.5") * (
                    vb[i] * wb[j] + wb[i] * vb[j]
                )
            else:
                raise ValueError(kind)
            mat[i, j] = mat[j, i] = weight * val
    return mat


def segments(rmax):
    pts = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    z = mp.mpf("0.4")
    while z < rmax:
        pts.append(z)
        z *= 2
    pts.append(rmax)
    return pts


def integrate(kind, c, s0, n, rmax, order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    out = mp.matrix(n)
    pts = segments(rmax)
    for a, b in zip(pts[:-1], pts[1:]):
        mid = mp.mpf("0.5") * (a + b)
        half = mp.mpf("0.5") * (b - a)
        for k in range(order):
            r = mid + half * nodes[k]
            out += half * weights[k] * contribution(kind, c, s0, r, n)
    return out, len(pts) - 1


def print_eigs(name, mat, limit):
    vals = mp.eigsy(mat, eigvals_only=True)
    print(f"  {name} eig low:")
    for val in vals[: min(limit, len(vals))]:
        print(f"    {mp.nstr(val, 25)}")
    return vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--kind", choices=["ka", "fgram", "kb"], default="kb")
    parser.add_argument("--rmax", type=str, default="12")
    parser.add_argument("--order", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    rmax = mp.mpf(args.rmax)
    mat, nseg = integrate(args.kind, c, s0, args.n, rmax, args.order)
    vals = print_eigs(args.kind, mat, args.n)
    print(
        f"endpoint confluent {args.kind} n={args.n} s0={s0} "
        f"dps={args.dps} rmax={rmax} order={args.order} "
        f"segments={nseg} min={mp.nstr(vals[0], 25)}"
    )


if __name__ == "__main__":
    main()
