#!/usr/bin/env python3
"""Feature-level identities for the endpoint-normalized P0 kernel.

Define

  v_lambda = (g_lambda-g_c)/(lambda-c),  lambda>c,
  g_lambda(x)=(1-lambda x)exp(-lambda(x-1)).

Then

  Ehat(lambda,mu)=<v_lambda,v_mu>,
  Phat0(lambda,mu)=1/2 log(lambda mu/c^2) Ehat(lambda,mu),
  d(lambda)=<v_lambda,g_c>,
  e(lambda)=zeta(lambda)/(lambda-c)=1/lambda.

This script verifies the closed-form normalized formulas against direct
endpoint quadrature.  It also prints endpoint/asymptotic row values to guide
the next analytic proof.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value
from endpoint_closed_form_p0 import beta_closed, e_kernel, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_feature_identity.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def v_value(c: float, lam: float, tau: float) -> float:
    return (g_value(lam, tau) - g_value(c, tau)) / (lam - c)


def direct_forms(args, lambdas):
    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts
    n = len(lambdas)
    ehat = np.zeros((n, n), dtype=float)
    drow = np.zeros(n, dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(args.c, float(tau))
        v = np.array([v_value(args.c, float(lam), float(tau)) for lam in lambdas])
        ehat += weight * np.outer(v, v)
        drow += weight * v * gc
    return sym(ehat), drow, segments, first


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--points", type=int, default=10)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=24)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    args = parser.parse_args()

    s_pts = np.linspace(args.T / args.points, args.T, args.points)
    lambdas = args.c * np.exp(s_pts)
    e_direct, d_direct, segments, first = direct_forms(args, lambdas)
    endpoint = lambdas - args.c
    e_closed = np.array(
        [
            [
                e_kernel(args.c, float(lam), float(mu)) / (dl * dm)
                for mu, dm in zip(lambdas, endpoint)
            ]
            for lam, dl in zip(lambdas, endpoint)
        ]
    )
    d_closed = np.array(
        [
            beta_closed(args.c, float(lam)) / (float(lam) - args.c)
            for lam in lambdas
        ]
    )
    phat = np.array(
        [
            [
                p0_kernel(args.c, float(lam), float(mu)) / (dl * dm)
                for mu, dm in zip(lambdas, endpoint)
            ]
            for lam, dl in zip(lambdas, endpoint)
        ]
    )
    print(
        f"endpoint normalized feature identity s=(0,{args.T:g}] points={args.points} "
        f"r=[0,{args.rmax:g}] segments={segments} first={first:.3e}"
    )
    print(f"  ||Ehat_closed-Ehat_direct|| = {np.linalg.norm(e_closed-e_direct):.12e}")
    print(f"  ||d_closed-d_direct||       = {np.linalg.norm(d_closed-d_direct):.12e}")
    vals = np.linalg.eigvalsh(sym(phat))
    print(
        f"  Phat point matrix min={vals[0]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -1e-10).sum():3d}"
    )
    print("  row samples:")
    for s, lam, d in zip(s_pts[: min(6, len(s_pts))], lambdas, d_closed):
        print(f"    s={s:.6g} e=1/lambda={1.0/lam:.12e} d={d:.12e}")


if __name__ == "__main__":
    main()
