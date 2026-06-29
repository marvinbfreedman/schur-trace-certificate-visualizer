#!/usr/bin/env python3
"""Range identity for the endpoint-normalized P0 theorem.

Let

  F_lambda(tau)=(exp(-lambda tau)-exp(-c tau))/(lambda-c),
  B=x d/dtau+1,  x=1+tau.

Then v_lambda=B F_lambda=(g_lambda-g_c)/(lambda-c), and

  Phat0(lambda,mu)=1/2(s_lambda+s_mu)<B F_lambda,B F_mu>.

For a finite combination F=sum alpha_i F_lambda_i and its spectral-log
companion H=sum alpha_i s_i F_lambda_i, the quadratic form is

  Phat0(alpha)=<B F,B H>.

The two rows are Hilbert pairings

  e(alpha)=<B F,h_z>,   h_z=c x^(-3/2)log x,
  d(alpha)=<B F,g_c>.

This script verifies the range identity directly.  It is the analytic form of
the remaining proof target: show <BF,BH> >= 0 on the range where these two
boundary pairings vanish.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value
from endpoint_closed_form_p0 import beta_closed, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_range_identity.py requires numpy") from exc


def f_value(c: float, lam: float, tau: float) -> float:
    return (math.exp(-lam * tau) - math.exp(-c * tau)) / (lam - c)


def bf_value(c: float, lam: float, tau: float) -> float:
    return (g_value(lam, tau) - g_value(c, tau)) / (lam - c)


def h_z(c: float, tau: float) -> float:
    x = 1.0 + tau
    return c * (x ** -1.5) * math.log(x)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=26)
    parser.add_argument("--ratio", type=float, default=1.35)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    s = np.linspace(args.T / args.n, args.T, args.n)
    lam = args.c * np.exp(s)
    alpha = rng.normal(size=args.n)
    # Give the vector enough sign variation to exercise the indefinite form.
    alpha -= alpha.mean()

    first = min(1e-4, 0.05 / float(lam[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    direct = 0.0
    e_row = 0.0
    d_row = 0.0
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        bf = np.array([bf_value(args.c, float(l), float(tau)) for l in lam])
        bh = s * bf
        bf_sum = float(alpha @ bf)
        bh_sum = float(alpha @ bh)
        direct += weight * bf_sum * bh_sum
        e_row += weight * bf_sum * h_z(args.c, float(tau))
        d_row += weight * bf_sum * g_value(args.c, float(tau))

    closed = 0.0
    for i, li in enumerate(lam):
        for j, lj in enumerate(lam):
            closed += (
                alpha[i]
                * alpha[j]
                * p0_kernel(args.c, float(li), float(lj))
                / ((li - args.c) * (lj - args.c))
            )
    e_closed = float(alpha @ (1.0 / lam))
    d_closed = float(
        alpha
        @ np.array(
            [beta_closed(args.c, float(l)) / (float(l) - args.c) for l in lam]
        )
    )

    print(
        f"endpoint normalized range identity n={args.n} s=(0,{args.T:g}] "
        f"r=[0,{args.rmax:g}] segments={segments}"
    )
    print(f"  <BF,BH> direct = {direct:.12e}")
    print(f"  Phat closed    = {closed:.12e}")
    print(f"  defect         = {direct-closed:.12e}")
    print(f"  e row direct   = {e_row:.12e}")
    print(f"  e row closed   = {e_closed:.12e}")
    print(f"  d row direct   = {d_row:.12e}")
    print(f"  d row closed   = {d_closed:.12e}")


if __name__ == "__main__":
    main()
