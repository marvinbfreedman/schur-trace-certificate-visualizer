#!/usr/bin/env python3
"""Generalized spectrum of Phat relative to the positive Ehat Gram.

For v_lambda=(g_lambda-g_c)/(lambda-c),

  Ehat(lambda,mu)=<v_lambda,v_mu>,
  Phat(lambda,mu)=1/2(s_lambda+s_mu)Ehat(lambda,mu).

This script studies the generalized eigenproblem Phat f = theta Ehat f and
the effect of imposing the two boundary rows

  e(lambda)=1/lambda=<v_lambda,h_z>,
  d(lambda)=<v_lambda,g_c>.

If the negative directions are endpoint modes, they should appear as a small
finite-dimensional negative generalized spectrum removed by e,d.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import beta_closed, e_kernel, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_generalized_spectrum.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def build(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s = 0.5 * args.T * (nodes + 1.0)
    w = 0.5 * args.T * weights
    root = np.sqrt(w)
    lam = args.c * np.exp(s)
    endpoint = lam - args.c
    ehat = np.zeros((args.order, args.order), dtype=float)
    phat = np.zeros_like(ehat)
    for i, x in enumerate(lam):
        for j, y in enumerate(lam[: i + 1]):
            denom = endpoint[i] * endpoint[j]
            eg = root[i] * e_kernel(args.c, float(x), float(y)) * root[j] / denom
            pg = root[i] * p0_kernel(args.c, float(x), float(y)) * root[j] / denom
            ehat[i, j] = ehat[j, i] = eg
            phat[i, j] = phat[j, i] = pg
    erow = root / lam
    drow = root * np.array(
        [beta_closed(args.c, float(x)) / (float(x) - args.c) for x in lam]
    )
    return s, lam, sym(ehat), sym(phat), erow, drow


def generalized(mat, gram, tol):
    vals_g, vecs_g = np.linalg.eigh(sym(gram))
    keep = vals_g > tol
    q = vecs_g[:, keep] / np.sqrt(vals_g[keep])
    reduced = sym(q.T @ sym(mat) @ q)
    vals, vecs = np.linalg.eigh(reduced)
    modes = q @ vecs
    return vals, modes, vals_g, keep


def restricted_basis(rows, n):
    constraints = np.vstack(rows)
    _, singular, vt = np.linalg.svd(constraints, full_matrices=True)
    rank = int((singular > 1e-12).sum())
    return vt[rank:].T, singular


def row_moments(mode, rows):
    return [float(row @ mode) for row in rows]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--show", type=int, default=6)
    args = parser.parse_args()

    s, lam, ehat, phat, erow, drow = build(args)
    vals, modes, gvals, keep = generalized(phat, ehat, args.tol)
    basis, singular = restricted_basis([erow, drow], args.order)
    rvals, _, rgvals, _ = generalized(basis.T @ phat @ basis, basis.T @ ehat @ basis, args.tol)
    print(
        f"endpoint normalized generalized spectrum lambda=[c,c exp({args.T:g})] "
        f"order={args.order}"
    )
    print(
        f"  Ehat gram min={gvals[0]:.12e} kept={keep.sum()} "
        f"dropped={(~keep).sum()}"
    )
    print(
        f"  full generalized min={vals[0]:.12e} second={vals[1]:.12e} "
        f"neg={(vals < -args.tol).sum():3d}"
    )
    print(
        f"  row-null generalized min={rvals[0]:.12e} "
        f"neg={(rvals < -args.tol).sum():3d}"
    )
    print("  row singular values: " + " ".join(f"{x:.12e}" for x in singular))
    print("  low generalized modes:")
    for k in range(min(args.show, len(vals))):
        moments = row_moments(modes[:, k], [erow, drow])
        # Convert weighted coefficients to a rough s-location diagnostic.
        coeff = modes[:, k]
        mass = coeff * coeff
        mass_sum = float(mass.sum())
        s_mean = float((mass @ s) / mass_sum) if mass_sum else float("nan")
        print(
            f"    {k}: theta={vals[k]: .12e} "
            f"<e>={moments[0]: .6e} <d>={moments[1]: .6e} "
            f"s_mean={s_mean:.6g}"
        )


if __name__ == "__main__":
    main()
