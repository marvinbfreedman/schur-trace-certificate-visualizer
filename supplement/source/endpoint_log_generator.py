#!/usr/bin/env python3
"""Scalar logarithmic-generator representation for H.

For boundary-null Laplace ranges

  F(tau)=sum_i a_i (exp(-lambda_i tau)-exp(-c tau)),
  H(tau)=sum_i a_i log(lambda_i/c) exp(-lambda_i tau),

the logarithm has the exact scalar representation

  H(tau)=int_0^inf [exp(-c u)F(tau)-F(tau+u)] du/u.

This makes C=E(F,H) a logarithmic generator form for the positive Sturm energy

  E(F,G)=int x^(7/2)F'G' - (3/2)int x^(3/2)FG.
"""

from __future__ import annotations

import argparse
import math
import random

from endpoint_boundary_null_kernel import composite_r_quadrature
from positive_branch_perturbation import quadrature


def values(c: float, lambdas: list[float], coeffs: list[float], tau: float):
    f = 0.0
    fp = 0.0
    h = 0.0
    hp = 0.0
    ec = math.exp(-c * tau)
    for lam, coeff in zip(lambdas, coeffs):
        el = math.exp(-lam * tau)
        ell = math.log(lam / c)
        f += coeff * (el - ec)
        fp += coeff * (-lam * el + c * ec)
        h += coeff * ell * el
        hp += coeff * ell * (-lam * el)
    return f, fp, h, hp


def shifted_log_value(c: float, lambdas: list[float], coeffs: list[float], tau: float, umax: float, order: int):
    pts, wts = quadrature(umax, order)
    total = 0.0
    for u, w in zip(pts, wts):
        uf = float(u)
        f_tau, _, _, _ = values(c, lambdas, coeffs, tau)
        f_shift, _, _, _ = values(c, lambdas, coeffs, tau + uf)
        total += float(w) * (math.exp(-c * uf) * f_tau - f_shift) / uf
    return total


def sturm_e(c: float, lambdas: list[float], coeffs: list[float], rmax: float, segment_order: int):
    r_pts, r_wts, _ = composite_r_quadrature(rmax, 1e-5, 1.35, segment_order)
    tau_pts = [math.expm1(float(r)) for r in r_pts]
    tau_wts = [math.exp(float(r)) * float(w) for r, w in zip(r_pts, r_wts)]
    total = 0.0
    for tau, w in zip(tau_pts, tau_wts):
        x = 1.0 + tau
        f, fp, h, hp = values(c, lambdas, coeffs, tau)
        total += w * (x ** 3.5 * fp * hp - 1.5 * x ** 1.5 * f * h)
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--umax", type=float, default=20.0)
    parser.add_argument("--u-order", type=int, default=600)
    parser.add_argument("--rmax", type=float, default=16.0)
    parser.add_argument("--segment-order", type=int, default=24)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    ss = sorted(rng.random() * args.T for _ in range(args.n))
    lambdas = [args.c * math.exp(s) for s in ss]
    coeffs = [rng.uniform(-1.0, 1.0) for _ in ss]

    print("endpoint logarithmic generator")
    print("  s nodes:", " ".join(f"{s:.6g}" for s in ss))
    print("  coeffs :", " ".join(f"{a:.6g}" for a in coeffs))
    max_defect = 0.0
    for tau in (0.0, 0.02, 0.1, 0.7, 2.0):
        _, _, h, _ = values(args.c, lambdas, coeffs, tau)
        approx = shifted_log_value(args.c, lambdas, coeffs, tau, args.umax, args.u_order)
        defect = approx - h
        max_defect = max(max_defect, abs(defect))
        print(
            f"  tau={tau:g}: H={h:.12e} shift={approx:.12e} "
            f"defect={defect:.12e}"
        )
    print(f"  max shift identity defect = {max_defect:.12e}")
    print(f"  C=E(F,H)                 = {sturm_e(args.c, lambdas, coeffs, args.rmax, args.segment_order):.12e}")


if __name__ == "__main__":
    main()
