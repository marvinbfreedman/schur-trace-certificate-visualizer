#!/usr/bin/env python3
"""Source-row form of the normalized endpoint constraints.

For B=x d/dtau+1 in L^2(x^(3/2)dtau),

  B^*h = -x h' - 3h/2.

The normalized rows become source moments of F:

  e(lambda)=<B F_lambda,h_z> = -c int F_lambda(tau) dtau = 1/lambda,

  d(lambda)=<B F_lambda,g_c>
           = <F_lambda, B^*g_c>,

where

  B^*g_c = (-c^2 x^2 + (7c/2)x - 3/2) exp(-c(x-1)).

This script verifies these source identities directly.  They are the two
ordinary moment constraints in the normalized Hardy theorem.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_closed_form_p0 import beta_closed

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_source_rows.py requires numpy") from exc


def f_value(c: float, lam: float, tau: float) -> float:
    return (math.exp(-lam * tau) - math.exp(-c * tau)) / (lam - c)


def bstar_gc(c: float, tau: float) -> float:
    x = 1.0 + tau
    return (-c * c * x * x + 3.5 * c * x - 1.5) * math.exp(-c * tau)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--points", type=int, default=8)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=28)
    parser.add_argument("--ratio", type=float, default=1.35)
    args = parser.parse_args()

    s_pts = np.linspace(args.T / args.points, args.T, args.points)
    lambdas = args.c * np.exp(s_pts)
    first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    max_e_defect = 0.0
    max_d_defect = 0.0
    print(
        f"endpoint normalized source rows s=(0,{args.T:g}] points={args.points} "
        f"r=[0,{args.rmax:g}] segments={segments}"
    )
    for s, lam in zip(s_pts, lambdas):
        e_int = 0.0
        d_int = 0.0
        for tau, tau_wt in zip(tau_pts, tau_wts):
            x = 1.0 + float(tau)
            f = f_value(args.c, float(lam), float(tau))
            e_int += float(tau_wt) * f
            d_int += float(tau_wt) * (x**1.5) * f * bstar_gc(args.c, float(tau))
        e_source = -args.c * e_int
        e_closed = 1.0 / float(lam)
        d_source = d_int
        d_closed = beta_closed(args.c, float(lam)) / (float(lam) - args.c)
        max_e_defect = max(max_e_defect, abs(e_source - e_closed))
        max_d_defect = max(max_d_defect, abs(d_source - d_closed))
        print(
            f"  s={s:.6g} e_src={e_source:.12e} e={e_closed:.12e} "
            f"d_src={d_source:.12e} d={d_closed:.12e}"
        )
    print(f"  max e defect = {max_e_defect:.12e}")
    print(f"  max d defect = {max_d_defect:.12e}")


if __name__ == "__main__":
    main()
