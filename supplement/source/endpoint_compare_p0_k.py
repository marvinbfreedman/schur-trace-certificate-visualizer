#!/usr/bin/env python3
"""Compare the finite-rank corrected P0 form with the boundary-null K form.

The previous endpoint branch produced the positive candidate

  K = L + C_A

with A=x d/dtau+3/2.  Since for boundary-null functions the A and B Dirichlet
Grams agree, it is natural to ask whether the new corrected anti-commutator

  P0 + M(beta^2+zeta^2)

is just K plus a positive remainder.  This script tests that possible shortcut.
If neither difference has a stable sign, the P0 Hardy inequality is genuinely
separate from the K positivity theorem.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature, feature
from endpoint_closed_form_certificate import build as build_closed
from endpoint_closed_form_certificate import sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_compare_p0_k.py requires numpy") from exc


def build_k(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    t_pts = 0.5 * args.T * (nodes + 1.0)
    t_wts = 0.5 * args.T * weights
    lambdas = args.c * np.exp(t_pts)
    root = np.sqrt(t_wts)

    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    kmat = np.zeros((args.order, args.order), dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        af = np.zeros(args.order, dtype=float)
        ah = np.zeros(args.order, dtype=float)
        for i, lam in enumerate(lambdas):
            af[i], ah[i], _, _ = feature(args.c, float(lam), float(tau))
        layer = math.log(x) * np.outer(af, af) + 0.5 * (
            np.outer(af, ah) + np.outer(ah, af)
        )
        kmat += weight * root[:, None] * layer * root[None, :]
    return sym(kmat), segments, first


def eig_line(name, mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<16} min={vals[0]: .12e} second={vals[1]: .12e} "
        f"max={vals[-1]: .12e} neg={(vals < -tol).sum():3d}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=17.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--M", type=float, default=14.0)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    p0, cmat, beta, zeta = build_closed(args)
    k, segments, first = build_k(args)
    corrected = p0 + args.M * (np.outer(beta, beta) + np.outer(zeta, zeta))

    print(
        f"endpoint P0/K comparison lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    eig_line("P0+rows", corrected, args.tol)
    eig_line("K", k, args.tol)
    eig_line("P0+rows-K", corrected - k, args.tol)
    eig_line("K-(P0+rows)", k - corrected, args.tol)
    eig_line("C+rows-K", cmat + args.M * (np.outer(beta, beta) + np.outer(zeta, zeta)) - k, args.tol)


if __name__ == "__main__":
    main()
