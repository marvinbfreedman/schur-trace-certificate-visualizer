#!/usr/bin/env python3
"""Mixed spectral derivative diagnostic for the boundary-null kernel.

Let

  K0(lambda,mu)=int x^(3/2) q_lambda q_mu
    log(x sqrt(lambda mu)/c) dx,

where q_lambda=(3/2-lambda x) exp(-lambda(x-1)).  The boundary-null kernel is
the conditionalization of K0 at lambda=c.  If

  J(lambda,mu)=partial_lambda partial_mu K0(lambda,mu)

is positive definite, then K is positive definite because

  K(lambda,mu)=int_c^lambda int_c^mu J(a,b) db da.

This script tests J directly with the endpoint-layer composite quadrature.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_mixed_derivative_kernel.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def q_and_dq(c: float, lambdas, x: float):
    tau = x - 1.0
    expv = np.exp(-lambdas * tau)
    q = (1.5 - lambdas * x) * expv
    dq = (lambdas * x * tau - x - 1.5 * tau) * expv
    return q, dq


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=16.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    nodes, _ = np.polynomial.legendre.leggauss(args.order)
    t_pts = 0.5 * args.T * (nodes + 1.0)
    lambdas = args.c * np.exp(t_pts)
    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    mat = np.zeros((args.order, args.order), dtype=float)
    for tau, w in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        logx = math.log(x)
        q, dq = q_and_dq(args.c, lambdas, x)
        ell = np.log(lambdas / args.c)
        # d_lambda d_mu of q_l q_m [log x + (ell_l+ell_m)/2].
        left = q / lambdas + ell * dq
        layer = (
            logx * np.outer(dq, dq)
            + 0.5 * np.outer(dq, left)
            + 0.5 * np.outer(left, dq)
        )
        mat += float(w) * (x ** 1.5) * layer

    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"endpoint mixed derivative J lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"segment_order={args.segment_order}"
    )
    print(
        f"  min={vals[0]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -args.tol).sum()}"
    )
    print("  low eigenvalues:", " ".join(f"{v:.3e}" for v in vals[:8]))


if __name__ == "__main__":
    main()
