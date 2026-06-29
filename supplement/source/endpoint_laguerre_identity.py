#!/usr/bin/env python3
"""Verify the endpoint Laguerre ODE and Green identity.

For

  h(u)=exp(5u/4)(c exp(u)-3/2) exp(-c exp(u))

write

  phi(u)=exp(3u/2) exp(-c exp(u)).

Then h(u)=-exp(-u/4) phi'(u).  For finite translates,

  A(r)=sum_i a_i h(t_i+r),
  B(r)=sum_i a_i (t_i+r) h(t_i+r),

put d_i=a_i exp(-t_i/4),

  Phi(r)=sum_i d_i phi(t_i+r),
  Psi(r)=sum_i d_i (t_i+r) phi(t_i+r).

The exact identity is

  A(r)B(r)
    = d/dr[-1/2 exp(-r/2) Phi(r)^2]
      + exp(-r/2) [Phi'(r) Psi'(r)-1/4 Phi(r)^2].

This is the closest local differential identity: the remaining term is a
constrained range term, not a pointwise square.
"""

from __future__ import annotations

import argparse
import math
import random

from positive_branch_perturbation import quadrature


def h_value(c: float, u: float) -> float:
    z = c * math.exp(u)
    return math.exp(1.25 * u) * (z - 1.5) * math.exp(-z)


def h_derivatives(c: float, u: float) -> tuple[float, float, float]:
    z = c * math.exp(u)
    h = h_value(c, u)
    # log h derivative for direct stable differentiation.
    m = 1.25 - z + z / (z - 1.5)
    mp = -z - 1.5 * z / ((z - 1.5) ** 2)
    return h, m * h, (mp + m * m) * h


def phi_values(c: float, u: float) -> tuple[float, float]:
    z = c * math.exp(u)
    phi = math.exp(1.5 * u) * math.exp(-z)
    return phi, (1.5 - z) * phi


def combo_values(c: float, ts: list[float], coeffs: list[float], r: float):
    A = 0.0
    B = 0.0
    Phi = 0.0
    Phip = 0.0
    Psi = 0.0
    Psip = 0.0
    for t, coeff in zip(ts, coeffs):
        u = t + r
        h = h_value(c, u)
        A += coeff * h
        B += coeff * u * h

        d = coeff * math.exp(-0.25 * t)
        phi, phip = phi_values(c, u)
        Phi += d * phi
        Phip += d * phip
        Psi += d * u * phi
        Psip += d * (phi + u * phip)
    return A, B, Phi, Phip, Psi, Psip


def green_density(c: float, ts: list[float], coeffs: list[float], r: float) -> float:
    _, _, Phi, Phip, _, Psip = combo_values(c, ts, coeffs, r)
    return math.exp(-0.5 * r) * (Phip * Psip - 0.25 * Phi * Phi)


def boundary_value(c: float, ts: list[float], coeffs: list[float], r: float) -> float:
    _, _, Phi, _, _, _ = combo_values(c, ts, coeffs, r)
    return -0.5 * math.exp(-0.5 * r) * Phi * Phi


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--order", type=int, default=500)
    parser.add_argument("--seed", type=int, default=8)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    ts = sorted(rng.random() * 3.0 for _ in range(5))
    coeffs = [rng.uniform(-1.0, 1.0) for _ in ts]

    print("endpoint Laguerre ODE / Green identity")
    print("  t nodes:", " ".join(f"{t:.6g}" for t in ts))
    print("  coeffs :", " ".join(f"{a:.6g}" for a in coeffs))

    max_ode = 0.0
    for k in range(41):
        u = args.rmax * k / 40.0
        z = args.c * math.exp(u)
        h, hp, hpp = h_derivatives(args.c, u)
        residual = hpp + (z - 2.0) * hp + (1.25 * z + 15.0 / 16.0) * h
        max_ode = max(max_ode, abs(residual))
    print(f"  max ODE residual on [0,{args.rmax:g}] = {max_ode:.12e}")

    pts, wts = quadrature(args.rmax, args.order)
    direct = 0.0
    residual = 0.0
    min_residual_density = float("inf")
    min_residual_at = None
    max_point_defect = 0.0
    for r, w in zip(pts, wts):
        rf = float(r)
        A, B, Phi, Phip, _, Psip = combo_values(args.c, ts, coeffs, rf)
        # Differentiate the boundary numerically only for the pointwise check.
        eps = 1e-5
        bd = (boundary_value(args.c, ts, coeffs, rf + eps) - boundary_value(args.c, ts, coeffs, rf - eps)) / (
            2.0 * eps
        )
        gd = math.exp(-0.5 * rf) * (Phip * Psip - 0.25 * Phi * Phi)
        if gd < min_residual_density:
            min_residual_density = gd
            min_residual_at = rf
        max_point_defect = max(max_point_defect, abs(A * B - bd - gd))
        direct += float(w) * A * B
        residual += float(w) * green_density(args.c, ts, coeffs, rf)

    b0 = boundary_value(args.c, ts, coeffs, 0.0)
    b1 = boundary_value(args.c, ts, coeffs, args.rmax)
    green_total = (b1 - b0) + residual
    print(f"  direct integral                 = {direct:.12e}")
    print(f"  boundary contribution           = {(b1 - b0):.12e}")
    print(f"  constrained residual integral   = {residual:.12e}")
    print(
        f"  min residual density            = {min_residual_density:.12e} "
        f"at r={min_residual_at:.6g}"
    )
    print(f"  Green identity total            = {green_total:.12e}")
    print(f"  integral defect                 = {direct - green_total:.12e}")
    print(f"  max pointwise defect            = {max_point_defect:.12e}")


if __name__ == "__main__":
    main()
