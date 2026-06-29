#!/usr/bin/env python3
"""Boundary-row representation for the P0 finite-rank correction.

Let x=1+tau and B=x d/dtau+1.  In the weighted space

  <f,g> = int_1^inf x^(3/2) f(x) g(x) dx,

the adjoint satisfies

  B^* h = -x h' - 3h/2.

The quotient rows for Lemma A are

  beta(lambda)=<u_lambda,g_c>,
  zeta(lambda)=1-c/lambda,
  u_lambda=B(e^(-lambda tau)-e^(-c tau)).

The second row is also a Hilbert-space pairing.  Since

  B^*(c x^(-3/2) log x) = -c x^(-3/2),

and F_lambda(0)=0,

  <u_lambda, c x^(-3/2) log x>
    = <B F_lambda, c x^(-3/2) log x>
    = -c int_0^inf F_lambda(tau) dtau
    = 1-c/lambda.

This converts the finite-rank target into

  P0(U) + 14(<U,g_c>^2 + <U,h_z>^2) >= 0,
  h_z(x)=c x^(-3/2) log x.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value


def h_zeta(c: float, tau: float) -> float:
    x = 1.0 + tau
    return c * x ** -1.5 * math.log(x)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--points", type=int, default=9)
    parser.add_argument("--rmax", type=float, default=20.0)
    parser.add_argument("--segment-order", type=int, default=40)
    parser.add_argument("--first-step", type=float, default=1e-6)
    parser.add_argument("--ratio", type=float, default=1.35)
    args = parser.parse_args()

    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, args.first_step, args.ratio, args.segment_order
    )
    tau_pts = [math.expm1(float(r)) for r in r_pts]
    tau_wts = [math.exp(float(r)) * float(w) for r, w in zip(r_pts, r_wts)]

    max_defect = 0.0
    print(
        f"endpoint P0 boundary rows s=[0,{args.T:g}] points={args.points} "
        f"r=[0,{args.rmax:g}] segments={segments}"
    )
    for i in range(1, args.points + 1):
        s = args.T * i / args.points
        lam = args.c * math.exp(s)
        zeta = 1.0 - args.c / lam
        pairing = 0.0
        beta = 0.0
        for tau, wt in zip(tau_pts, tau_wts):
            x = 1.0 + tau
            u = g_value(lam, tau) - g_value(args.c, tau)
            gc = g_value(args.c, tau)
            weight = wt * x**1.5
            pairing += weight * u * h_zeta(args.c, tau)
            beta += weight * u * gc
        defect = pairing - zeta
        max_defect = max(max_defect, abs(defect))
        print(
            f"  s={s:.6g} beta={beta:.12e} "
            f"zeta={zeta:.12e} <u,h_z>={pairing:.12e} "
            f"defect={defect:.3e}"
        )
    print(f"  max zeta pairing defect = {max_defect:.12e}")


if __name__ == "__main__":
    main()
