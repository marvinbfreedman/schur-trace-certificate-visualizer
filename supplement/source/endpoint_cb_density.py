#!/usr/bin/env python3
"""Scan the paired Cauchy-Binet density for the endpoint kernel P.

With

  f_0(t,r)=h(t+r),    f_1(t,r)=(t+r)h(t+r),

the kernel

  P(t,s)=1/2 int [f_1(t,r)f_0(s,r)+f_0(t,r)f_1(s,r)] dr

has Cauchy-Binet density

  S(T,S;R)=sum_alpha det[f_{1-alpha_j}(t_i,r_j)]
                    det[f_{alpha_j}(s_i,r_j)].

For principal minors T=S this becomes the paired density

  S(T;R)=sum_alpha det[f_{1-alpha_j}(t_i,r_j)]
                  det[f_{alpha_j}(t_i,r_j)].

If S(T;R)>=0 pointwise for ordered T,R, then the principal minors of P are
nonnegative by Cauchy-Binet even though individual mixed type-word signs fail.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random


def h_value(c: float, u: float) -> float:
    eu = math.exp(u)
    return math.exp(1.25 * u) * (c * eu - 1.5) * math.exp(-c * eu)


def f_value(c: float, t: float, r: float, kind: int) -> float:
    u = t + r
    value = h_value(c, u)
    if kind == 0:
        return value
    if kind == 1:
        return u * value
    raise ValueError(kind)


def det(matrix: list[list[float]]) -> float:
    n = len(matrix)
    a = [row[:] for row in matrix]
    sign = 1.0
    out = 1.0
    for i in range(n):
        pivot = max(range(i, n), key=lambda k: abs(a[k][i]))
        if abs(a[pivot][i]) == 0.0:
            return 0.0
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            sign *= -1.0
        piv = a[i][i]
        out *= piv
        for k in range(i + 1, n):
            factor = a[k][i] / piv
            for j in range(i + 1, n):
                a[k][j] -= factor * a[i][j]
    return sign * out


def typed_det(c: float, rows: tuple[float, ...], cols: tuple[float, ...], kinds) -> float:
    return det(
        [
            [f_value(c, float(t), float(r), kinds[j]) for j, r in enumerate(cols)]
            for t in rows
        ]
    )


def paired_density(c: float, rows: tuple[float, ...], cols: tuple[float, ...]) -> float:
    n = len(rows)
    total = 0.0
    for kinds in itertools.product((0, 1), repeat=n):
        complement = tuple(1 - k for k in kinds)
        total += typed_det(c, rows, cols, complement) * typed_det(c, rows, cols, kinds)
    return total


def linspace(a: float, b: float, n: int) -> list[float]:
    if n == 1:
        return [a]
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def random_tuple(rng: random.Random, n: int, scale: float) -> tuple[float, ...]:
    return tuple(sorted(rng.random() * scale for _ in range(n)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=5.0)
    parser.add_argument("--grid", type=int, default=10)
    parser.add_argument("--max-order", type=int, default=7)
    parser.add_argument("--sample-limit", type=int, default=2000)
    parser.add_argument("--random", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--tol", type=float, default=1e-18)
    args = parser.parse_args()

    grid = linspace(0.0, args.T, args.grid)
    rng = random.Random(args.seed)
    print(
        f"endpoint paired Cauchy-Binet density T=[0,{args.T:g}] "
        f"grid={args.grid} max_order={args.max_order}"
    )
    for n in range(1, args.max_order + 1):
        worst = float("inf")
        wrong = 0
        total = 0
        worst_rows = None
        worst_cols = None

        if args.random:
            iterator = (
                (random_tuple(rng, n, args.T), random_tuple(rng, n, args.T))
                for _ in range(args.random)
            )
        else:
            iterator = itertools.product(
                itertools.combinations(grid, n), itertools.combinations(grid, n)
            )

        for rows, cols in iterator:
            value = paired_density(args.c, rows, cols)
            if value < worst:
                worst = value
                worst_rows = rows
                worst_cols = cols
            wrong += value < -args.tol
            total += 1
            if not args.random and total >= args.sample_limit:
                break

        print(
            f"  order={n}: min_density={worst:.12e} wrong={wrong}/{total} "
            f"rows={worst_rows} cols={worst_cols}"
        )


if __name__ == "__main__":
    main()
