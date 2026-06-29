#!/usr/bin/env python3
r"""High-frequency obstruction for the compact commuted-kernel model.

The finite certificate used

    S_m^comm(f) = sum_{r=0}^m <K D^r f, D^r f>,

with K the endpoint-B integral kernel.  If K is treated as a compact
continuum operator on L^2(0,L), this form cannot dominate the full W^{m,2}
norm on an infinite-dimensional space.  A smooth packet supported away from
the moving endpoint trace interval is in ker R, has W^{m,2} norm one, and
converges weakly at top derivative order, forcing the compact quadratic form
to zero.

This script gives a direct numerical version of that obstruction using
compactly supported polynomial-sine packets.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import endpoint_b_quadrature, endpoint_b_vw  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def poly_mul(a, b):
    out = [mp.mpf("0") for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def shifted_power(root, power):
    return [mp.binomial(power, j) * ((-root) ** (power - j)) for j in range(power + 1)]


def poly_derivative(poly, deriv):
    if deriv >= len(poly):
        return [mp.mpf("0")]
    out = []
    for k in range(deriv, len(poly)):
        out.append(poly[k] * mp.factorial(k) / mp.factorial(k - deriv))
    return out


def poly_value(poly, x):
    total = mp.mpf("0")
    for coeff in reversed(poly):
        total = total * x + coeff
    return total


def sin_derivative(freq, phase, deriv):
    amp = freq**deriv
    mod = deriv % 4
    if mod == 0:
        return amp * mp.sin(phase)
    if mod == 1:
        return amp * mp.cos(phase)
    if mod == 2:
        return -amp * mp.sin(phase)
    return -amp * mp.cos(phase)


def packet_derivative(poly, a, b, freq, deriv, x):
    if x <= a or x >= b:
        return mp.mpf("0")
    phase = freq * (x - a)
    total = mp.mpf("0")
    for j in range(deriv + 1):
        dp = poly_derivative(poly, j)
        total += (
            mp.binomial(deriv, j)
            * poly_value(dp, x)
            * sin_derivative(freq, phase, deriv - j)
        )
    return total


def legendre_interval(a, b, order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    mid = (a + b) / 2
    half = (b - a) / 2
    return [mid + half * x for x in nodes], [half * w for w in weights]


def packet_profile(args, freq):
    a = mp.mpf(args.support_a)
    b = mp.mpf(args.support_b)
    power = args.vanish_order
    left = shifted_power(a, power)
    right = shifted_power(b, power)
    # (s-a)^p (b-s)^p = (s-a)^p (-1)^p (s-b)^p.
    poly = poly_mul(left, [((-1) ** power) * c for c in right])
    s_nodes, s_weights = legendre_interval(a, b, args.s_order)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.rmax), args.r_order)

    deriv_values = []
    wnorm = mp.mpf("0")
    for deriv in range(args.order + 1):
        vals = [packet_derivative(poly, a, b, freq, deriv, s) for s in s_nodes]
        deriv_values.append(vals)
        wnorm += mp.fsum(w * val * val for w, val in zip(s_weights, vals))

    energy = mp.mpf("0")
    for deriv, vals in enumerate(deriv_values):
        kquad = mp.mpf("0")
        for r, rw in zip(r_nodes, r_weights):
            vint = mp.mpf("0")
            wint = mp.mpf("0")
            for s, sw, val in zip(s_nodes, s_weights, vals):
                v, wfun = endpoint_b_vw(s, r, mp.pi)
                vint += sw * v * val
                wint += sw * wfun * val
            kquad += rw * mp.e ** (mp.mpf("2.5") * r) * (
                r * vint * vint + vint * wint
            )
        energy += kquad
    return energy, wnorm, energy / wnorm if wnorm else mp.mpf("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=10)
    parser.add_argument("--frequencies", default="8 12 16 20 24 32")
    parser.add_argument("--support-a", default="0.8")
    parser.add_argument("--support-b", default="1.8")
    parser.add_argument("--vanish-order", type=int, default=12)
    parser.add_argument("--s-order", type=int, default=180)
    parser.add_argument("--r-order", type=int, default=18)
    parser.add_argument("--rmax", default="12")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="commuted_compact_obstruction.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print("Compact commuted-kernel obstruction")
    print("  freq       S_m          W_m        S_m/W_m")
    for piece in args.frequencies.replace(",", " ").split():
        freq = mp.mpf(piece)
        energy, wnorm, ratio = packet_profile(args, freq)
        rows.append(
            {
                "frequency": f(freq),
                "energy": f(energy),
                "sobolevNorm2": f(wnorm),
                "ratio": f(ratio),
            }
        )
        print(
            f"  {fmt(freq, 5):>5} {fmt(energy, 8):>12} "
            f"{fmt(wnorm, 8):>12} {fmt(ratio, 8):>12}",
            flush=True,
        )

    data = {
        "order": args.order,
        "support": [f(mp.mpf(args.support_a)), f(mp.mpf(args.support_b))],
        "vanishOrder": args.vanish_order,
        "rows": rows,
        "interpretation": (
            "For the compact kernel model S_m^comm, smooth oscillatory packets "
            "supported away from the endpoint trace interval drive S_m/W_m "
            "toward zero. The continuum elliptic estimate therefore needs the "
            "local positive terms from the true commuted Sturm identity."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
