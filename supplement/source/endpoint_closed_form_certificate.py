#!/usr/bin/env python3
"""Closed-form finite-rank certificates for P0 and C.

This uses the closed formulas from endpoint_closed_form_p0.py to build
Legendre-Galerkin sections in s=log(lambda/c), without tau quadrature:

  P0(lambda,mu)=1/2(s_lambda+s_mu)E(lambda,mu),
  C=P0+R2,

where E, beta, and zeta are all evaluated by the J_p recurrences.  It then
tests the explicit finite-rank lower bound

  P0 + M(beta^2+zeta^2) >= 0,
  C  + M(beta^2+zeta^2) >= 0.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import (
    beta_closed,
    c_kernel,
    p0_kernel,
    zeta_closed,
)

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_closed_form_certificate.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def corrected_values(mat, beta, zeta, m):
    correction = np.outer(beta, beta) + np.outer(zeta, zeta)
    vals = np.linalg.eigvalsh(sym(mat + m * correction))
    return vals


def threshold_m(mat, beta, zeta, tol):
    lo = 0.0
    hi = 1.0
    while corrected_values(mat, beta, zeta, hi)[0] < -tol:
        hi *= 2.0
        if hi > 1e9:
            return float("inf")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if corrected_values(mat, beta, zeta, mid)[0] < -tol:
            lo = mid
        else:
            hi = mid
    return hi


def build(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)

    p0 = np.zeros((args.order, args.order), dtype=float)
    cmat = np.zeros_like(p0)
    for i, lam in enumerate(lambdas):
        for j, mu in enumerate(lambdas[: i + 1]):
            p = root[i] * p0_kernel(args.c, float(lam), float(mu)) * root[j]
            q = root[i] * c_kernel(args.c, float(lam), float(mu)) * root[j]
            p0[i, j] = p0[j, i] = p
            cmat[i, j] = cmat[j, i] = q
    beta = root * np.array([beta_closed(args.c, float(lam)) for lam in lambdas])
    zeta = root * np.array([zeta_closed(args.c, float(lam)) for lam in lambdas])
    return sym(p0), sym(cmat), beta, zeta


def report(name, mat, beta, zeta, args):
    vals = np.linalg.eigvalsh(sym(mat))
    corrected = corrected_values(mat, beta, zeta, args.M)
    mstar = threshold_m(mat, beta, zeta, args.tol)
    print(
        f"  {name}: min={vals[0]:.12e} second={vals[1]:.12e} "
        f"neg={(vals < -args.tol).sum():3d}"
    )
    print(
        f"    threshold M ~= {mstar:.12e}; "
        f"with M={args.M:g}: min={corrected[0]:.12e} "
        f"neg={(corrected < -args.tol).sum():3d}"
    )
    print("    low corrected: " + " ".join(f"{v:.3e}" for v in corrected[:8]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--M", type=float, default=14.0)
    args = parser.parse_args()

    p0, cmat, beta, zeta = build(args)
    print(
        f"endpoint closed-form certificate lambda=[c,c exp({args.T:g})] "
        f"order={args.order}"
    )
    report("P0", p0, beta, zeta, args)
    report("C ", cmat, beta, zeta, args)


if __name__ == "__main__":
    main()
