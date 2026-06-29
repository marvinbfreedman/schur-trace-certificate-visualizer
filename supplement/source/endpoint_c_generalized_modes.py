#!/usr/bin/env python3
"""Generalized negative modes of the endpoint cross form C.

The completed Sturm identity gives

  C(F,H)=int_0^inf x^(3/2) ((xF)')((xH)') dtau,

and the positive reference Dirichlet form is

  E(F,G)=int_0^inf x^(3/2) ((xF)')((xG)') dtau.

On the boundary-null Laplace range

  F_lambda(tau)=exp(-lambda tau)-exp(-c tau),
  H_lambda(tau)=log(lambda/c) exp(-lambda tau),

the two-negative-square problem for C should appear as a finite negative
generalized spectrum of

  C v = rho E v.

This script extracts those E-relative modes.  It is a numerical model for the
desired analytic decomposition C = P - D_2, where P >= 0 and rank(D_2) <= 2.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "endpoint_c_generalized_modes.py requires numpy; run with python"
    ) from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def f_pair(c: float, lam: float, tau: float):
    """Return F, F', H, H' for the boundary-null feature."""
    ell = math.log(lam / c)
    ec = math.exp(-c * tau)
    el = math.exp(-lam * tau)
    f = el - ec
    fp = -lam * el + c * ec
    h = ell * el
    hp = -lam * ell * el
    return f, fp, h, hp


def build_forms(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)

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
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        ux_f = np.zeros(args.order, dtype=float)
        ux_h = np.zeros(args.order, dtype=float)
        for i, lam in enumerate(lambdas):
            f, fp, h, hp = f_pair(args.c, float(lam), float(tau))
            ux_f[i] = f + x * fp
            ux_h[i] = h + x * hp
        weight = float(tau_wt) * x**1.5
        emat += weight * root[:, None] * np.outer(ux_f, ux_f) * root[None, :]
        cmat += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(ux_f, ux_h) + np.outer(ux_h, ux_f))
            * root[None, :]
        )

    return s_pts, s_wts, lambdas, root, sym(emat), sym(cmat), segments, first


def generalized_spectrum(emat, cmat, keep_rel):
    evals, evecs = np.linalg.eigh(sym(emat))
    cutoff = keep_rel * max(1.0, float(evals[-1]))
    keep = evals > cutoff
    if not keep.any():
        raise SystemExit("E whitening kept no modes; lower --keep-rel")
    white = evecs[:, keep] / np.sqrt(evals[keep])[None, :]
    generator = sym(white.T @ sym(cmat) @ white)
    rhos, modes = np.linalg.eigh(generator)
    original_modes = white @ modes
    return evals, keep, cutoff, rhos, original_modes, generator


def physical_coeff(mode, root, s_wts):
    coeff = mode / root
    norm = math.sqrt(max(float(np.sum(s_wts * coeff * coeff)), 0.0))
    if norm > 0.0:
        coeff = coeff / norm
    return coeff


def mode_report(c, s_pts, s_wts, lambdas, root, rhos, modes, tol, show):
    shown = 0
    rows = []
    for j, rho in enumerate(rhos):
        if rho >= -tol and shown >= show:
            break
        if rho >= -tol:
            continue
        coeff = physical_coeff(modes[:, j], root, s_wts)
        moments = [float(np.sum(s_wts * (s_pts**k) * coeff)) for k in range(6)]
        lam_moments = [
            float(np.sum(s_wts * ((lambdas - c) ** k) * coeff))
            for k in range(1, 4)
        ]
        # Endpoint Taylor coefficients of F(tau)=int coeff(s)F_lambda(tau) ds.
        fp0 = float(np.sum(s_wts * coeff * (c - lambdas)))
        fpp0 = float(np.sum(s_wts * coeff * (lambdas * lambdas - c**2)))
        rows.append(
            "    mode {idx}: rho={rho:.12e} "
            "m0={m0: .6e} m1={m1: .6e} m2={m2: .6e} "
            "m3={m3: .6e} fp0={fp0: .6e} fpp0={fpp0: .6e} "
            "lam1={l1: .6e} lam2={l2: .6e} lam3={l3: .6e}".format(
                idx=j,
                rho=float(rho),
                m0=moments[0],
                m1=moments[1],
                m2=moments[2],
                m3=moments[3],
                fp0=fp0,
                fpp0=fpp0,
                l1=lam_moments[0],
                l2=lam_moments[1],
                l3=lam_moments[2],
            )
        )
        shown += 1
        if shown >= show:
            break
    return rows


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
    parser.add_argument("--keep-rel", type=float, default=1e-14)
    parser.add_argument("--show-modes", type=int, default=4)
    args = parser.parse_args()

    s_pts, s_wts, lambdas, root, emat, cmat, segments, first = build_forms(args)
    evals, keep, cutoff, rhos, modes, _ = generalized_spectrum(
        emat, cmat, args.keep_rel
    )
    cvals = np.linalg.eigvalsh(sym(cmat))
    neg_rhos = int((rhos < -args.tol).sum())

    print(
        f"endpoint C generalized modes lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    print(
        f"  E min={evals[0]:.12e} max={evals[-1]:.12e} "
        f"kept={int(keep.sum())}/{len(evals)} cutoff={cutoff:.3e}"
    )
    print(
        f"  ordinary C min={cvals[0]:.12e} max={cvals[-1]:.12e} "
        f"neg={(cvals < -args.tol).sum()}"
    )
    print(
        f"  generalized rho min={rhos[0]:.12e} max={rhos[-1]:.12e} "
        f"neg={neg_rhos}"
    )
    low = " ".join(f"{x:.6e}" for x in rhos[: min(8, len(rhos))])
    print(f"  low generalized rhos: {low}")
    if args.show_modes:
        print("  E-relative negative mode diagnostics:")
        for line in mode_report(
            args.c,
            s_pts,
            s_wts,
            lambdas,
            root,
            rhos,
            modes,
            args.tol,
            args.show_modes,
        ):
            print(line)


if __name__ == "__main__":
    main()
