#!/usr/bin/env python3
"""Formal Taylor-jet tests for the finite-core reduced Volterra/Weyl kernel.

This is the full finite-core object

  K_red(s,t) = integral_0^inf [u+(s+t)/2] cosh(omega[u+(s+t)/2])
                 A_s(u) A_t(u) du,

where

  A_s(u) = Psi(s+u) / Psi(s),
  Psi(v) = tilde Phi_3(v/2)

or the corresponding raw finite core.  The script builds the confluent matrix

  [ d_s^i d_t^j K_red(s0,t0) / (i! j!) ]_{i,j=0}^{n-1},  t0=s0,

directly from formal Taylor series at the integrand level.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from mp_partial_k_corr import weights_for  # noqa: E402


def zero(n):
    return [mp.mpf("0") for _ in range(n)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [k * x for x in a]


def mul(a, b):
    n = min(len(a), len(b))
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


def cosh_series(a):
    ep = exp_series(a)
    em = exp_series(scale(a, -1))
    return scale(add(ep, em), mp.mpf("0.5"))


def pieces_for(kind):
    if kind == "raw1":
        mode_weights = {1: mp.mpf(1)}
    else:
        mode_weights = weights_for(kind)
    pieces = []
    for mode, weight in mode_weights.items():
        n2 = mode * mode
        c = mp.pi * n2
        pieces.append((weight * 4 * c * c, mp.mpf("2.25"), c))
        pieces.append((weight * -6 * c, mp.mpf("1.25"), c))
    return pieces


def psi_series(v0, pieces, n):
    ev0 = mp.e**v0
    e_delta = [1 / mp.factorial(k) for k in range(n)]
    out = zero(n)
    for coeff, beta, c in pieces:
        exponent = scale(e_delta, -c * ev0)
        exponent[0] += beta * v0
        if n > 1:
            exponent[1] += beta
        term = scale(exp_series(exponent), coeff)
        out = add(out, term)
    return out


def weight_series(center, omega, n):
    # Series in z = delta_s + delta_t for
    # (center+z/2) cosh(omega(center+z/2)).
    r = zero(n)
    r[0] = center
    if n > 1:
        r[1] = mp.mpf("0.5")
    return mul(r, cosh_series(scale(r, omega)))


def ratio_series(s0, u, pieces, denom, n):
    return div(psi_series(s0 + u, pieces, n), denom)


def contribution(s0, u, omega, pieces, denom, n):
    a = ratio_series(s0, u, pieces, denom, n)
    w = weight_series(s0 + u, omega, 2 * n - 1)
    mat = mp.matrix(n)
    for i in range(n):
        for j in range(i + 1):
            val = mp.mpf("0")
            for p in range(i + 1):
                for q in range(j + 1):
                    k = (i - p) + (j - q)
                    val += a[p] * a[q] * w[k] * mp.binomial(k, i - p)
            mat[i, j] = mat[j, i] = val
    return mat


def segments(umax):
    pts = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    z = mp.mpf("0.4")
    while z < umax:
        pts.append(z)
        z *= 2
    pts.append(umax)
    return pts


def integrate(s0, omega, pieces, n, umax, order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    denom = psi_series(s0, pieces, n)
    out = mp.matrix(n)
    pts = segments(umax)
    for a, b in zip(pts[:-1], pts[1:]):
        mid = mp.mpf("0.5") * (a + b)
        half = mp.mpf("0.5") * (b - a)
        for k in range(order):
            u = mid + half * nodes[k]
            out += half * weights[k] * contribution(s0, u, omega, pieces, denom, n)
    return out, len(pts) - 1, denom[0]


def print_eigs(mat, limit):
    vals = mp.eigsy(mat, eigvals_only=True)
    print("  K_red jet eig low:")
    for val in vals[: min(limit, len(vals))]:
        print(f"    {mp.nstr(val, 25)}")
    return vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--n", type=int, default=9)
    parser.add_argument("--s0", default="0.5")
    parser.add_argument("--umax", default="10")
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
    omega = mp.mpf(args.omega)
    umax = mp.mpf(args.umax)
    pieces = pieces_for(args.kind)
    mat, nseg, psi0 = integrate(s0, omega, pieces, args.n, umax, args.order)
    vals = print_eigs(mat, args.n)
    print(
        f"finite-core Volterra confluent kind={args.kind} omega={omega} "
        f"n={args.n} s0={s0} dps={args.dps} umax={umax} "
        f"order={args.order} segments={nseg}"
    )
    print(f"  Psi(s0)={mp.nstr(psi0, 25)}")
    print(f"  min={mp.nstr(vals[0], 25)}")


if __name__ == "__main__":
    main()
