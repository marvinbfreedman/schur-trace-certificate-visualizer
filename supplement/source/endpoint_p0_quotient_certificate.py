#!/usr/bin/env python3
"""Codimension-two quotient certificate candidate for Lemma A.

Lemma A asks for

  ind_-(P0) <= 2,
  P0(lambda,mu)=1/2(s_lambda+s_mu)<u_lambda,u_mu>,
  s_lambda=log(lambda/c).

Numerically the correct quotient is not given by the raw moments (1,s).  The
stable codimension-two candidate is

  beta(lambda) = <u_lambda,g_c>,
  zeta(lambda) = 1 - c/lambda.

The second functional has a simple Laplace meaning:

  zeta(lambda) = -c int_0^inf F_lambda(tau) dtau.

This script verifies that P0 restricted to

  beta(alpha)=0,  zeta(alpha)=0

is nonnegative up to the quadrature/nullspace floor.  Proving this quotient
positivity analytically proves Lemma A.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "endpoint_p0_quotient_certificate.py requires numpy; run with python"
    ) from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def restricted_eigs(mat, rows, tol):
    constraints = np.vstack(rows)
    _, singular, vt = np.linalg.svd(constraints, full_matrices=True)
    rank = int((singular > 1e-10).sum())
    basis = vt[rank:].T
    restricted = sym(basis.T @ sym(mat) @ basis)
    vals = np.linalg.eigvalsh(restricted)
    return vals, singular


def build(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    logs = np.log(lambdas / args.c)
    zeta = 1.0 - args.c / lambdas

    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    p0 = np.zeros((args.order, args.order), dtype=float)
    beta = np.zeros(args.order, dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(args.c, float(tau))
        u = np.array(
            [g_value(float(lam), float(tau)) - gc for lam in lambdas],
            dtype=float,
        )
        gram = weight * root[:, None] * np.outer(u, u) * root[None, :]
        p0 += 0.5 * (logs[:, None] + logs[None, :]) * gram
        beta += weight * root * u * gc

    return s_pts, s_wts, root, lambdas, sym(p0), beta, root * zeta, segments, first


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

    _, _, _, _, p0, beta, zeta, segments, first = build(args)
    vals = np.linalg.eigvalsh(sym(p0))
    qvals, singular = restricted_eigs(p0, [beta, zeta], args.tol)
    print(
        f"endpoint P0 quotient certificate lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    print(
        f"  P0 min={vals[0]:.12e} second={vals[1]:.12e} "
        f"neg={(vals < -args.tol).sum():3d}"
    )
    print(
        "  quotient beta=<u,g_c>, zeta=1-c/lambda: "
        f"min={qvals[0]:.12e} neg={(qvals < -args.tol).sum():3d} "
        f"neg1e-10={(qvals < -1e-10).sum():3d} dim={len(qvals)}"
    )
    print(
        "  constraint singular values: "
        + " ".join(f"{x:.6e}" for x in singular[: min(4, len(singular))])
    )
    print("  low quotient eigenvalues: " + " ".join(f"{v:.3e}" for v in qvals[:8]))


if __name__ == "__main__":
    main()
