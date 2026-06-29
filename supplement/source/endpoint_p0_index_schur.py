#!/usr/bin/env python3
"""Index transfer from P0 to C=P0+R2.

The exact split

  C = P0 + R2,
  P0(lambda,mu)=1/2(log(lambda/c)+log(mu/c)) E(lambda,mu),
  rank(R2)<=2,

suggests the index theorem:

  ind_-(P0) <= 2, and R2 does not increase the index.

This script tests the second clause in the sharp form needed for a proof.  Let
N0 be the negative spectral subspace of P0.  Then C has at most two negative
squares if C is nonnegative on N0^perp, with the usual Moore-Penrose range
condition for the nullspace of that block.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_p0_index_schur.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def build_split(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    ells = np.log(lambdas / args.c)

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
    bvec = np.zeros(args.order, dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(args.c, float(tau))
        u = np.array(
            [g_value(float(lam), float(tau)) - gc for lam in lambdas],
            dtype=float,
        )
        ell_g = np.array(
            [ell * g_value(float(lam), float(tau)) for ell, lam in zip(ells, lambdas)],
            dtype=float,
        )
        emat += weight * root[:, None] * np.outer(u, u) * root[None, :]
        cmat += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(u, ell_g) + np.outer(ell_g, u))
            * root[None, :]
        )
        bvec += weight * root * u * gc

    p0 = 0.5 * (ells[:, None] + ells[None, :]) * emat
    ell_weighted = root * ells
    r2 = 0.5 * (np.outer(bvec, ell_weighted) + np.outer(ell_weighted, bvec))
    return s_pts, s_wts, sym(emat), sym(p0), sym(r2), sym(cmat), segments, first


def block_over_negative(total, reference, tol):
    vals, vecs = np.linalg.eigh(sym(reference))
    neg = vals < -tol
    nspace = vecs[:, neg]
    comp = vecs[:, ~neg]
    total = sym(total)
    aa = nspace.T @ total @ nspace
    ab = nspace.T @ total @ comp
    bb = comp.T @ total @ comp
    bb_vals, bb_vecs = np.linalg.eigh(sym(bb))
    keep = bb_vals > tol
    null = ~keep
    range_defect = 0.0
    if null.any():
        range_defect = float(np.linalg.norm(ab @ bb_vecs[:, null]))
    if keep.any():
        inv = (bb_vecs[:, keep] / bb_vals[keep]) @ bb_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
    else:
        schur = aa
    return vals, bb_vals, np.linalg.eigvalsh(sym(schur)), range_defect


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

    _, _, emat, p0, r2, cmat, segments, first = build_split(args)
    print(
        f"endpoint P0 index transfer lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    for name, mat in (("E", emat), ("P0", p0), ("R2", r2), ("C", cmat)):
        vals = np.linalg.eigvalsh(sym(mat))
        print(
            f"  {name:<2} min={vals[0]: .12e} max={vals[-1]: .12e} "
            f"neg={(vals < -args.tol).sum():3d}"
        )

    p0_vals, c_bb, c_schur, c_defect = block_over_negative(cmat, p0, args.tol)
    print(
        "  C on P0-nonnegative block: "
        f"min={c_bb[0]:.12e} neg={(c_bb < -args.tol).sum():3d} "
        f"null={(abs(c_bb) <= args.tol).sum():3d}"
    )
    print(
        "  Schur C over neg(P0): "
        f"min={c_schur[0]:.12e} max={c_schur[-1]:.12e} "
        f"range_defect={c_defect:.12e}"
    )
    print(
        "  low P0 eigenvalues: "
        + " ".join(f"{x:.3e}" for x in p0_vals[: min(8, len(p0_vals))])
    )


if __name__ == "__main__":
    main()
