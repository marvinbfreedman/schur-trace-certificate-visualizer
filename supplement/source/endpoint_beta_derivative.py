#!/usr/bin/env python3
"""Matrix diagnostics for the beta-derivative endpoint Hardy form.

Let

  Phi_beta(r)=sum_i a_i exp(beta(t_i+r)) exp(-c exp(t_i+r)).

At beta=3/2, Psi=partial_beta Phi_beta.  The Green identity from
endpoint_laguerre_identity.py gives

  Q = 1/2 Phi(0)^2
      + 1/2 partial_beta int e^{-r/2} |Phi_beta'(r)|^2 dr
      - 1/4 int e^{-r/2} |Phi(r)|^2 dr.

Equivalently,

  2Q = Boundary + D_beta - 1/2 N.

This script builds Galerkin matrices for those three pieces on translate nodes
t_i and reports which piece carries the endpoint positivity.
"""

from __future__ import annotations

import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_beta_derivative.py requires numpy; run with python") from exc


def phi_beta(c: float, beta: float, u: float) -> float:
    return math.exp(beta * u - c * math.exp(u))


def phi_beta_prime(c: float, beta: float, u: float) -> float:
    return (beta - c * math.exp(u)) * phi_beta(c, beta, u)


def d_beta_phi_prime(c: float, beta: float, u: float) -> float:
    # partial_beta [partial_u phi_beta] = phi_beta + u partial_u phi_beta.
    phi = phi_beta(c, beta, u)
    return phi + u * (beta - c * math.exp(u)) * phi


def sym(mat):
    return 0.5 * (mat + mat.T)


def report(name: str, mat, tol: float):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<18} min={vals[0]: .12e} max={vals[-1]: .12e} "
        f"neg={(vals < -tol).sum():3d}"
    )
    return vals


def build(args):
    beta = args.beta
    t_pts, t_wts = quadrature(args.T, args.t_order)
    r_pts, r_wts = quadrature(args.rmax, args.r_order)
    root = np.sqrt(t_wts)

    nmat = np.zeros((args.t_order, args.t_order), dtype=float)
    dmat = np.zeros_like(nmat)
    dbeta = np.zeros_like(nmat)
    boundary = np.zeros_like(nmat)

    bvals = np.array([phi_beta(args.c, beta, float(t)) for t in t_pts])
    boundary = root[:, None] * np.outer(bvals, bvals) * root[None, :]

    for r, w in zip(r_pts, r_wts):
        rf = float(r)
        weight = math.exp(-0.5 * rf)
        phi = np.array([phi_beta(args.c, beta, float(t + r)) for t in t_pts])
        phip = np.array([phi_beta_prime(args.c, beta, float(t + r)) for t in t_pts])
        dphip = np.array([d_beta_phi_prime(args.c, beta, float(t + r)) for t in t_pts])
        nmat += float(w) * weight * root[:, None] * np.outer(phi, phi) * root[None, :]
        dmat += float(w) * weight * root[:, None] * np.outer(phip, phip) * root[None, :]
        dbeta += (
            float(w)
            * weight
            * root[:, None]
            * (np.outer(dphip, phip) + np.outer(phip, dphip))
            * root[None, :]
        )

    q2 = boundary + dbeta - 0.5 * nmat
    q = 0.5 * q2
    free_hardy = dmat + 0.25 * boundary - (1.0 / 16.0) * nmat
    return {
        "N": nmat,
        "D": dmat,
        "D_beta": dbeta,
        "Boundary": boundary,
        "D_beta-.5N": dbeta - 0.5 * nmat,
        "2Q": q2,
        "Q": q,
        "free_Hardy": free_hardy,
        "t_pts": t_pts,
        "t_wts": t_wts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--beta", type=float, default=1.5)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--t-order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--r-order", type=int, default=300)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    mats = build(args)
    print(
        f"endpoint beta-derivative split beta={args.beta:g} c={args.c:g} "
        f"t=[0,{args.T:g}] order={args.t_order} r=[0,{args.rmax:g}]"
    )
    for name in ("N", "D", "D_beta", "Boundary", "D_beta-.5N", "2Q", "Q", "free_Hardy"):
        report(name, mats[name], args.tol)

    # Check how much endpoint boundary is needed on the negative subspace of
    # D_beta-.5N.
    vals, vecs = np.linalg.eigh(sym(mats["D_beta-.5N"]))
    mask = vals < -args.tol
    if mask.any():
        block = vecs[:, mask].T @ sym(mats["Boundary"]) @ vecs[:, mask]
        scale = np.diag(1.0 / np.sqrt(-vals[mask]))
        ratios = np.linalg.eigvalsh(sym(scale @ block @ scale))
        print(
            "  Boundary/-(D_beta-.5N) on negative subspace: "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e} dim={mask.sum()}"
        )

    bvals, bvecs = np.linalg.eigh(sym(mats["Boundary"]))
    bmask = bvals <= args.tol
    if bmask.any():
        null_basis = bvecs[:, bmask]
        restricted = null_basis.T @ sym(mats["D_beta-.5N"]) @ null_basis
        rvals = np.linalg.eigvalsh(sym(restricted))
        print(
            "  (D_beta-.5N) on Boundary-null: "
            f"min={rvals[0]:.12e} max={rvals[-1]:.12e} "
            f"dim={len(rvals)}"
        )


if __name__ == "__main__":
    main()
