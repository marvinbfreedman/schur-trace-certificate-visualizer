#!/usr/bin/env python3
"""Finite-rank lower-bound form for the cross-kernel C.

The same rows that certify the P0 quotient also certify the C quotient:

  beta(lambda)=<u_lambda,g_c>,
  zeta(lambda)=1-c/lambda=<u_lambda,h_z>.

Numerically,

  C + M(beta beta^T + zeta zeta^T) >= 0

with threshold M about 11.78 on the reliable spectral windows.  The rounded
M=14 bound is compatible with the P0 finite-rank target and gives the same
two-negative-square conclusion for C.
"""

from __future__ import annotations

import argparse
import math

from endpoint_p0_commutator_identity import build_forms
from endpoint_p0_finite_rank_bound import min_eig_for_m, threshold_m
from endpoint_p0_quotient_certificate import build, sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_c_finite_rank_bound.py requires numpy") from exc


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
    parser.add_argument("--M", type=float, default=14.0)
    args = parser.parse_args()

    _, _, _, _, _, beta, zeta, segments, first = build(args)
    cmat = build_forms(args)["C"]
    cvals = np.linalg.eigvalsh(sym(cmat))
    mstar = threshold_m(cmat, beta, zeta, args.tol)
    low_m, vals_m = min_eig_for_m(cmat, beta, zeta, args.M)
    print(
        f"endpoint C finite-rank bound lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    print(
        f"  C min={cvals[0]:.12e} second={cvals[1]:.12e} "
        f"neg={(cvals < -args.tol).sum():3d}"
    )
    print(f"  scalar threshold M ~= {mstar:.12e}")
    print(
        f"  with M={args.M:g}: min={low_m:.12e} "
        f"neg={(vals_m < -args.tol).sum():3d}"
    )
    print("  low corrected eigenvalues: " + " ".join(f"{v:.3e}" for v in vals_m[:8]))


if __name__ == "__main__":
    main()
