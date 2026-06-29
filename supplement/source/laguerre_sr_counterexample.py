#!/usr/bin/env python3
"""Exact Wronskian obstruction to the proposed Laguerre SR theorem.

For K(X,Y)=(XY-a)exp(-XY), the x-derivatives at fixed X give

  d_X^k K(X,Y) = (-1)^k Y^k (XY-a-k) exp(-XY).

After coalescing the Y columns, the local n-by-n sign is controlled by the
Wronskian of q_k(z)=z^k(z-a-k), k=0,...,n-1.  The proposed reverse-regular
signature would require this Wronskian to be positive throughout the domain.

At a=3/2, n=7, z=16/5 (> pi), the Wronskian is exactly negative, so the
confluent Laguerre sign-regularity theorem cannot hold even on the endpoint
rectangle X>=1, Y>=pi.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations


def trim(poly: list[Fraction]) -> list[Fraction]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    n = max(len(p), len(q))
    return trim(
        [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]
    )


def mul(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return trim(out)


def derivative(poly: list[Fraction], order: int) -> list[Fraction]:
    out = poly[:]
    for _ in range(order):
        out = [Fraction(i + 1) * out[i + 1] for i in range(len(out) - 1)] or [
            Fraction(0)
        ]
    return trim(out)


def det_poly(matrix: list[list[list[Fraction]]]) -> list[Fraction]:
    n = len(matrix)
    out = [Fraction(0)]
    for perm in permutations(range(n)):
        inversions = sum(
            1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j]
        )
        term = [Fraction(-1 if inversions % 2 else 1)]
        for i, j in enumerate(perm):
            term = mul(term, matrix[i][j])
        out = add(out, term)
    return trim(out)


def eval_poly(poly: list[Fraction], z: Fraction) -> Fraction:
    out = Fraction(0)
    for coeff in reversed(poly):
        out = out * z + coeff
    return out


def laguerre_wronskian(n: int, a: Fraction) -> list[Fraction]:
    polys = []
    for k in range(n):
        qk = [Fraction(0)] * (k + 2)
        qk[k] = -(a + k)
        qk[k + 1] = 1
        polys.append(qk)

    matrix = [[derivative(polys[j], i) for j in range(n)] for i in range(n)]
    return det_poly(matrix)


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=7)
    parser.add_argument("--a", type=parse_fraction, default=Fraction(3, 2))
    parser.add_argument("--z", type=parse_fraction, default=Fraction(16, 5))
    args = parser.parse_args()

    poly = laguerre_wronskian(args.n, args.a)
    value = eval_poly(poly, args.z)
    N = args.n * (args.n - 1) // 2

    print(f"n = {args.n}")
    print(f"a = {args.a}")
    print(f"z = {args.z}")
    print(f"signature epsilon_n = (-1)^{N} = {(-1) ** N}")
    print("Wronskian polynomial coefficients, increasing powers:")
    print("[" + ", ".join(str(c) for c in poly) + "]")
    print(f"W_{args.n - 1}({args.z}) = {value}")
    print(f"float(W) = {float(value):.12g}")
    if value < 0:
        print("Conclusion: local reverse sign-regularity fails.")
    else:
        print("Conclusion: this point does not obstruct local reverse sign-regularity.")


if __name__ == "__main__":
    main()
