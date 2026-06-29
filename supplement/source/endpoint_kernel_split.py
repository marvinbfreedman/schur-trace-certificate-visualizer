#!/usr/bin/env python3
"""Split the boundary-null kernel K into positive and cross pieces.

K = L + C, where

  L = int x^(3/2) log(x) (AF_lambda)(AF_mu) dx,
  C = 1/2 int x^(3/2) [
        (AF_lambda)(AH_mu) + (AF_mu)(AH_lambda)
      ] dx.

This checks whether the cross term has a finite negative index and whether the
positive L piece dominates its negative subspace.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature, feature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_kernel_split.py requires numpy; run with python") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def generalized_ratios(pos, ref, tol):
    vals, vecs = np.linalg.eigh(sym(ref))
    mask = vals < -tol
    if not mask.any():
        return vals, None
    block = vecs[:, mask].T @ sym(pos) @ vecs[:, mask]
    scale = np.diag(1.0 / np.sqrt(-vals[mask]))
    return vals, np.linalg.eigvalsh(sym(scale @ block @ scale))


def schur_relative(total, reference, tol):
    vals, vecs = np.linalg.eigh(sym(reference))
    mask = vals < -tol
    if not mask.any():
        return None, 0.0, 0
    neg = vecs[:, mask]
    comp = vecs[:, ~mask]
    total = sym(total)
    aa = neg.T @ total @ neg
    ab = neg.T @ total @ comp
    bb = comp.T @ total @ comp
    bb_vals, bb_vecs = np.linalg.eigh(sym(bb))
    keep = bb_vals > tol
    defect = 0.0
    if keep.any():
        inv = (bb_vecs[:, keep] / bb_vals[keep]) @ bb_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
    else:
        schur = aa
    if (~keep).any():
        defect = float(np.linalg.norm(ab @ bb_vecs[:, ~keep]))
    return schur, defect, int(keep.sum())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=16.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--show-modes", type=int, default=0)
    args = parser.parse_args()

    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    t_pts = 0.5 * args.T * (nodes + 1.0)
    t_wts = 0.5 * args.T * weights
    lambdas = args.c * np.exp(t_pts)
    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts
    root = np.sqrt(t_wts)

    lmat = np.zeros((args.order, args.order), dtype=float)
    cmat = np.zeros_like(lmat)
    for tau, w in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(w) * x ** 1.5
        af = np.zeros(args.order, dtype=float)
        ah = np.zeros(args.order, dtype=float)
        for i, lam in enumerate(lambdas):
            af[i], ah[i], _, _ = feature(args.c, float(lam), float(tau))
        lmat += weight * math.log(x) * root[:, None] * np.outer(af, af) * root[None, :]
        cmat += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(af, ah) + np.outer(ah, af))
            * root[None, :]
        )

    total = lmat + cmat
    print(
        f"endpoint K split lambda=[c,c exp({args.T:g})] order={args.order} "
        f"r=[0,{args.rmax:g}] segments={segments}"
    )
    for name, mat in (("L", lmat), ("C", cmat), ("K", total)):
        vals = np.linalg.eigvalsh(sym(mat))
        print(
            f"  {name:<2} min={vals[0]: .12e} max={vals[-1]: .12e} "
            f"neg={(vals < -args.tol).sum():3d}"
        )
    _, ratios = generalized_ratios(lmat, cmat, args.tol)
    if ratios is not None:
        print(
            "  L/(-C) on neg(C): "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e} dim={len(ratios)}"
        )
    schur, defect, rank = schur_relative(total, cmat, args.tol)
    if schur is not None:
        svals = np.linalg.eigvalsh(sym(schur))
        print(
            "  Schur K over neg(C): "
            f"min={svals[0]:.12e} max={svals[-1]:.12e} "
            f"bb_rank={rank} range_defect={defect:.12e}"
        )
    if args.show_modes:
        vals, vecs = np.linalg.eigh(sym(cmat))
        print("  negative C mode moments in s=log(lambda/c):")
        for j, val in enumerate(vals[: args.show_modes]):
            if val >= -args.tol:
                break
            coeff = vecs[:, j] / root
            norm = math.sqrt(float(np.sum(t_wts * coeff * coeff)))
            coeff = coeff / norm
            moments = [float(np.sum(t_wts * (t_pts ** k) * coeff)) for k in range(5)]
            print(
                f"    mode {j}: eig={val:.12e} "
                + " ".join(f"m{k}={moments[k]: .6e}" for k in range(5))
            )


if __name__ == "__main__":
    main()
