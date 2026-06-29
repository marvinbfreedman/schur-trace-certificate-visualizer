#!/usr/bin/env python3
"""Sturm-Liouville form of the cross-kernel C.

For boundary-null spectral features

  F_lambda(tau)=exp(-lambda tau)-exp(-c tau),
  H_lambda(tau)=log(lambda/c) exp(-lambda tau),
  A=3/2+x d/dtau, x=1+tau,

the cross form is

  C(F,H)=int x^(3/2) (AF)(AH) dtau.

Since F(0)=0, integration by parts gives the exact equivalent form

  C(F,H)=int x^(7/2) F' H' dtau - (3/2) int x^(3/2) F H dtau.

Equivalently, completing the first-order expression with

  U=d/dtau+1/x

gives

  C(F,H)=int x^(7/2) (UF)(UH) dtau
        =int x^(3/2) ((xF)')((xH)') dtau.

This script verifies the identity and scans simple endpoint/Taylor moment
restrictions for the negative-index problem.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature, feature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_c_sturm.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def null_min(mat, rows, tol):
    if not rows:
        vals = np.linalg.eigvalsh(sym(mat))
        return vals[0], int((vals < -tol).sum()), len(vals)
    constraints = np.vstack(rows)
    _, singular, vt = np.linalg.svd(constraints, full_matrices=True)
    rank = int((singular > 1e-10).sum())
    basis = vt[rank:].T
    restricted = basis.T @ sym(mat) @ basis
    vals = np.linalg.eigvalsh(sym(restricted))
    return vals[0], int((vals < -tol).sum()), len(vals)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=90)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    direct = np.zeros((args.order, args.order), dtype=float)
    sturm = np.zeros_like(direct)
    first_order = np.zeros_like(direct)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        af = np.zeros(args.order, dtype=float)
        ah = np.zeros(args.order, dtype=float)
        fprime = np.zeros(args.order, dtype=float)
        hprime = np.zeros(args.order, dtype=float)
        fval = np.zeros(args.order, dtype=float)
        hval = np.zeros(args.order, dtype=float)
        ux_f = np.zeros(args.order, dtype=float)
        ux_h = np.zeros(args.order, dtype=float)
        ec = math.exp(-args.c * float(tau))
        for i, lam in enumerate(lambdas):
            ell = math.log(float(lam) / args.c)
            el = math.exp(-float(lam) * float(tau))
            af[i], ah[i], _, _ = feature(args.c, float(lam), float(tau))
            fval[i] = el - ec
            hval[i] = ell * el
            fprime[i] = -float(lam) * el + args.c * ec
            hprime[i] = -float(lam) * ell * el
            ux_f[i] = fval[i] + x * fprime[i]
            ux_h[i] = hval[i] + x * hprime[i]
        direct += (
            0.5
            * float(tau_wt)
            * x**1.5
            * root[:, None]
            * (np.outer(af, ah) + np.outer(ah, af))
            * root[None, :]
        )
        sturm += (
            float(tau_wt)
            * root[:, None]
            * (
                x**3.5 * 0.5 * (np.outer(fprime, hprime) + np.outer(hprime, fprime))
                - 1.5 * x**1.5 * 0.5 * (np.outer(fval, hval) + np.outer(hval, fval))
            )
            * root[None, :]
        )
        first_order += (
            float(tau_wt)
            * x**1.5
            * root[:, None]
            * 0.5
            * (np.outer(ux_f, ux_h) + np.outer(ux_h, ux_f))
            * root[None, :]
        )

    defect = float(np.linalg.norm(direct - sturm))
    first_defect = float(np.linalg.norm(direct - first_order))
    vals = np.linalg.eigvalsh(sym(direct))
    print(
        f"endpoint C Sturm form lambda=[c,c exp({args.T:g})] order={args.order} "
        f"r=[0,{args.rmax:g}] segments={segments}"
    )
    print(f"  identity defect ||direct-sturm|| = {defect:.12e}")
    print(f"  first-order defect              = {first_defect:.12e}")
    print(
        f"  C min={vals[0]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -args.tol).sum()}"
    )

    candidates = {
        "1": np.ones_like(s_pts),
        "s": s_pts,
        "s^2": s_pts * s_pts,
        "lambda-c": lambdas - args.c,
        "s(lambda-c)": s_pts * (lambdas - args.c),
        "s(3/2-lambda)": s_pts * (1.5 - lambdas),
    }
    for names in (("1", "s"), ("1", "s^2"), ("s", "lambda-c"), ("s", "s(3/2-lambda)")):
        rows = [root * candidates[name] for name in names]
        low, neg, dim = null_min(direct, rows, args.tol)
        print(f"  null {names}: min={low:.12e} neg={neg} dim={dim}")


if __name__ == "__main__":
    main()
