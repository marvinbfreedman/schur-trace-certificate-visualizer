#!/usr/bin/env python3
"""Per-shift contribution in the logarithmic-generator representation.

The scalar formula

  H(tau)=int_0^inf [exp(-c u)F(tau)-F(tau+u)] du/u

gives shift contributions

  C_u(F)=exp(-c u)E(F,F)-E(F,F(.+u)).

This script shows that proving positivity or a two-negative-square bound for
each fixed u is too strong; the final theorem must use cancellation over u.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_shift_contribution.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def f_values(c: float, lam: float, tau: float):
    ec = math.exp(-c * tau)
    el = math.exp(-lam * tau)
    return el - ec, -lam * el + c * ec


def shift_matrix(c: float, lambdas, root, u: float, rmax: float, segment_order: int):
    first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, _ = composite_r_quadrature(rmax, first, 1.35, segment_order)
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts
    n = len(lambdas)
    energy = np.zeros((n, n), dtype=float)
    shifted = np.zeros_like(energy)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        f = np.zeros(n)
        fp = np.zeros(n)
        g = np.zeros(n)
        gp = np.zeros(n)
        for i, lam in enumerate(lambdas):
            f[i], fp[i] = f_values(c, float(lam), float(tau))
            g[i], gp[i] = f_values(c, float(lam), float(tau) + u)
        energy += (
            float(tau_wt)
            * root[:, None]
            * (x**3.5 * np.outer(fp, fp) - 1.5 * x**1.5 * np.outer(f, f))
            * root[None, :]
        )
        shifted += (
            float(tau_wt)
            * root[:, None]
            * (x**3.5 * np.outer(fp, gp) - 1.5 * x**1.5 * np.outer(f, g))
            * root[None, :]
        )
    return math.exp(-c * u) * energy - sym(shifted)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=70)
    parser.add_argument("--rmax", type=float, default=16.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)

    print(
        f"endpoint shift contribution lambda=[c,c exp({args.T:g})] "
        f"order={args.order}"
    )
    for u in (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0):
        mat = shift_matrix(args.c, lambdas, root, u, args.rmax, args.segment_order)
        vals = np.linalg.eigvalsh(sym(mat))
        print(
            f"  u={u:<5g} min={vals[0]: .12e} max={vals[-1]: .12e} "
            f"neg={(vals < -args.tol).sum():3d}"
        )


if __name__ == "__main__":
    main()
