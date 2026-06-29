#!/usr/bin/env python3
"""Schur closure of K=L+C using the beta,zeta rows.

After the quotient theorem

  C >= 0 on {beta=0,zeta=0},

the full positivity K=L+C can be checked by the Moore-Penrose Schur complement
relative to the same row space.  This script builds

  K = C + L,
  L = int log(1+tau) u_lambda u_mu x^(3/2) dtau,

and reports:

  1. K on the beta,zeta nullspace;
  2. the range defect;
  3. the Schur block over span{beta,zeta}.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value
from endpoint_p0_commutator_identity import build_forms
from endpoint_p0_quotient_certificate import build, restricted_eigs, sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_k_beta_zeta_schur.py requires numpy") from exc


def build_l(args, lambdas, root):
    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, _ = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts
    lmat = np.zeros((args.order, args.order), dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5 * math.log(x)
        gc = g_value(args.c, float(tau))
        u = np.array(
            [g_value(float(lam), float(tau)) - gc for lam in lambdas],
            dtype=float,
        )
        lmat += weight * root[:, None] * np.outer(u, u) * root[None, :]
    return sym(lmat)


def schur_over_rows(mat, rows, tol):
    constraints = np.vstack(rows)
    q, _ = np.linalg.qr(constraints.T)
    row_dim = constraints.shape[0]
    row_basis = q[:, :row_dim]
    _, singular, vt = np.linalg.svd(constraints, full_matrices=True)
    rank = int((singular > 1e-10).sum())
    null_basis = vt[rank:].T

    mat = sym(mat)
    aa = row_basis.T @ mat @ row_basis
    ab = row_basis.T @ mat @ null_basis
    bb = null_basis.T @ mat @ null_basis
    bb_vals, bb_vecs = np.linalg.eigh(sym(bb))
    keep = bb_vals > tol
    if keep.any():
        inv = (bb_vecs[:, keep] / bb_vals[keep]) @ bb_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
    else:
        schur = aa
    defect = 0.0
    if (~keep).any():
        defect = float(np.linalg.norm(ab @ bb_vecs[:, ~keep]))
    return bb_vals, np.linalg.eigvalsh(sym(schur)), defect, singular


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=100)
    parser.add_argument("--rmax", type=float, default=17.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    _, _, root, lambdas, _, beta, zeta, segments, first = build(args)
    cmat = build_forms(args)["C"]
    lmat = build_l(args, lambdas, root)
    kmat = sym(cmat + lmat)
    kvals = np.linalg.eigvalsh(kmat)
    qvals = restricted_eigs(kmat, [beta, zeta], args.tol)[0]
    bb_vals, schur_vals, defect, singular = schur_over_rows(
        kmat, [beta, zeta], args.tol
    )
    print(
        f"endpoint K beta/zeta Schur lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    print(
        f"  K min={kvals[0]:.12e} max={kvals[-1]:.12e} "
        f"neg={(kvals < -args.tol).sum():3d}"
    )
    print(
        f"  K on beta,zeta nullspace: min={qvals[0]:.12e} "
        f"neg={(qvals < -args.tol).sum():3d}"
    )
    print(
        f"  Schur over beta,zeta rows: min={schur_vals[0]:.12e} "
        f"max={schur_vals[-1]:.12e} range_defect={defect:.12e}"
    )
    print(
        "  constraint singular values: "
        + " ".join(f"{x:.6e}" for x in singular[: min(4, len(singular))])
    )
    print("  low quotient eigenvalues: " + " ".join(f"{v:.3e}" for v in qvals[:8]))


if __name__ == "__main__":
    main()
