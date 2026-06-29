#!/usr/bin/env python3
"""Exact rank-two source split for the endpoint cross form C.

Let

  g_lambda(tau)=((1+tau)e^(-lambda tau))'
               =(1-lambda(1+tau))e^(-lambda tau),
  u_lambda=g_lambda-g_c,
  ell_lambda=log(lambda/c).

Then the completed Sturm form gives

  E(lambda,mu)=<u_lambda,u_mu>,
  C(lambda,mu)=1/2(<u_lambda,ell_mu g_mu>
                 +<ell_lambda g_lambda,u_mu>).

Since g_lambda=u_lambda+g_c,

  C = P0 + R2,
  P0(lambda,mu)=1/2(ell_lambda+ell_mu)E(lambda,mu),
  R2(lambda,mu)=1/2 ell_mu <u_lambda,g_c>
               +1/2 ell_lambda <g_c,u_mu>.

The source R2 has rank at most two.  This split is useful, but P0 is not
positive by itself; the final proof still needs a commutator/Schur argument.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_c_rank2_split.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def g_value(lam: float, tau: float):
    x = 1.0 + tau
    return (1.0 - lam * x) * math.exp(-lam * tau)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=90)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    ells = np.log(lambdas / args.c)

    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    emat = np.zeros((args.order, args.order), dtype=float)
    cmat = np.zeros_like(emat)
    bvec = np.zeros(args.order, dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(args.c, float(tau))
        u = np.array(
            [g_value(float(lam), float(tau)) - gc for lam in lambdas],
            dtype=float,
        )
        ell_g = np.array(
            [ell * g_value(float(lam), float(tau)) for ell, lam in zip(ells, lambdas)],
            dtype=float,
        )
        emat += weight * root[:, None] * np.outer(u, u) * root[None, :]
        cmat += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(u, ell_g) + np.outer(ell_g, u))
            * root[None, :]
        )
        bvec += weight * root * u * gc

    p0 = 0.5 * (ells[:, None] + ells[None, :]) * emat
    # Weighted-coordinate form of
    # 1/2 ell_mu <u_lambda,g_c> + 1/2 ell_lambda <g_c,u_mu>.
    ell_weighted = root * ells
    r2 = 0.5 * (np.outer(bvec, ell_weighted) + np.outer(ell_weighted, bvec))
    defect = cmat - p0 - r2

    print(
        f"endpoint C rank-two split lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments}"
    )
    print(f"  split defect ||C-P0-R2|| = {np.linalg.norm(defect):.12e}")
    for name, mat in (("E", emat), ("P0", p0), ("R2", r2), ("C", cmat)):
        vals = np.linalg.eigvalsh(sym(mat))
        rank = int(np.linalg.matrix_rank(sym(mat), tol=args.tol))
        print(
            f"  {name:<2} min={vals[0]: .12e} max={vals[-1]: .12e} "
            f"neg={(vals < -args.tol).sum():3d} rank={rank:3d}"
        )
    rvals = np.linalg.eigvalsh(sym(r2))
    print(
        "  nonzero R2 eigenvalues: "
        + " ".join(f"{v:.12e}" for v in rvals[np.abs(rvals) > args.tol])
    )


if __name__ == "__main__":
    main()
