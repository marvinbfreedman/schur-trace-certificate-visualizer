#!/usr/bin/env python3
"""Boundary-null spectral kernel for the endpoint beta-Hardy form.

Normalize the endpoint range by

  lambda = c exp(t),      tau = exp(r)-1.

If Phi(0)=0, then after setting alpha_i=d_i phi(t_i) the range function is

  Phi(r) = x^(3/2) F(tau),  x=1+tau,
  F(tau) = sum_i alpha_i exp(-lambda_i tau),  sum_i alpha_i=0.

Using the base point lambda=c, write the boundary-null span with difference
features

  F_lambda(tau)=exp(-lambda tau)-exp(-c tau),
  H_lambda(tau)=log(lambda/c) exp(-lambda tau),

where H is the spectral-log companion.  With A=3/2+x d/dtau, the
boundary-null form is

  R/2 = int_0^inf x^(3/2) [
          log(x) (AF)^2 + (AF)(AH)
        ] dtau.

This script builds the symmetric kernel for that quadratic form and tests its
positivity on lambda in [c, c exp(T)].
"""

from __future__ import annotations

import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_boundary_null_kernel.py requires numpy; run with python") from exc


def feature(c: float, lam: float, tau: float) -> tuple[float, float, float, float]:
    x = 1.0 + tau
    e_lam = math.exp(-lam * tau)
    e_c = math.exp(-c * tau)
    log_lam = math.log(lam / c)

    f = e_lam - e_c
    fp = -lam * e_lam + c * e_c
    h = log_lam * e_lam
    hp = -lam * log_lam * e_lam

    af = 1.5 * f + x * fp
    ah = 1.5 * h + x * hp
    return af, ah, f, h


def sym(mat):
    return 0.5 * (mat + mat.T)


def composite_r_quadrature(end: float, first: float, ratio: float, order: int):
    edges = [0.0]
    x = first
    while x < end:
        edges.append(x)
        x *= ratio
    if edges[-1] < end:
        edges.append(end)

    all_pts = []
    all_wts = []
    nodes, weights = np.polynomial.legendre.leggauss(order)
    for left, right in zip(edges[:-1], edges[1:]):
        mid = 0.5 * (left + right)
        half = 0.5 * (right - left)
        all_pts.extend(mid + half * nodes)
        all_wts.extend(half * weights)
    return np.array(all_pts), np.array(all_wts), len(edges) - 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--taumax", type=float, default=200.0)
    parser.add_argument("--tau-order", type=int, default=500)
    parser.add_argument("--quad", choices=("tau", "r-composite"), default="r-composite")
    parser.add_argument("--rmax", type=float, default=16.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    t_pts, t_wts = quadrature(args.T, args.order)
    lambdas = args.c * np.exp(t_pts)
    if args.quad == "tau":
        tau_pts, tau_wts = quadrature(args.taumax, args.tau_order)
        quad_label = f"tau=[0,{args.taumax:g}] tau_order={args.tau_order}"
    else:
        first = args.first_step
        if first <= 0.0:
            first = min(1e-4, 0.05 / float(lambdas[-1]))
        r_pts, r_wts, segments = composite_r_quadrature(
            args.rmax, first, args.ratio, args.segment_order
        )
        tau_pts = np.expm1(r_pts)
        tau_wts = np.exp(r_pts) * r_wts
        quad_label = (
            f"r-composite=[0,{args.rmax:g}] segments={segments} "
            f"segment_order={args.segment_order} first={first:.3e}"
        )

    root = np.sqrt(t_wts)
    kmat = np.zeros((args.order, args.order), dtype=float)
    af_min = float("inf")
    for tau, w in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = (x ** 1.5) * float(w)
        af = np.zeros(args.order, dtype=float)
        ah = np.zeros(args.order, dtype=float)
        for i, lam in enumerate(lambdas):
            afi, ahi, _, _ = feature(args.c, float(lam), float(tau))
            af[i] = afi
            ah[i] = ahi
            af_min = min(af_min, afi)
        layer = math.log(x) * np.outer(af, af) + 0.5 * (
            np.outer(af, ah) + np.outer(ah, af)
        )
        kmat += weight * root[:, None] * layer * root[None, :]

    vals = np.linalg.eigvalsh(sym(kmat))
    print(
        f"endpoint boundary-null spectral kernel lambda=[c,c exp({args.T:g})] "
        f"order={args.order} {quad_label}"
    )
    print(
        f"  min={vals[0]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -args.tol).sum()}"
    )
    print("  low eigenvalues:", " ".join(f"{v:.3e}" for v in vals[:8]))


if __name__ == "__main__":
    main()
