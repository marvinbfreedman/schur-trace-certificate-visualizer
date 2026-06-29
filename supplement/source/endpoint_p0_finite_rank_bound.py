#!/usr/bin/env python3
"""Finite-rank lower-bound form of Lemma A.

The quotient certificate

  P0 >= 0 on {beta=0, zeta=0}

is equivalent to saying that P0 has a positive semidefinite completion after
adding a form in the two rows beta and zeta.  Numerically the simple scalar
completion

  P0 + M (beta beta^T + zeta zeta^T) >= 0

already works with M about 13.04 on the reliable spectral windows.  A clean
analytic target is therefore the slightly rounded inequality

  P0 + 14 (beta^2 + zeta^2) >= 0.

Proving this explicit finite-rank lower bound proves ind_-(P0) <= 2.
"""

from __future__ import annotations

import argparse
import math

from endpoint_p0_quotient_certificate import build, sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_p0_finite_rank_bound.py requires numpy") from exc


def min_eig_for_m(p0, beta, zeta, m):
    rows = np.outer(beta, beta) + np.outer(zeta, zeta)
    vals = np.linalg.eigvalsh(sym(p0 + m * rows))
    return vals[0], vals


def threshold_m(p0, beta, zeta, tol):
    lo = 0.0
    hi = 1.0
    while min_eig_for_m(p0, beta, zeta, hi)[0] < -tol:
        hi *= 2.0
        if hi > 1e8:
            return float("inf")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if min_eig_for_m(p0, beta, zeta, mid)[0] < -tol:
            lo = mid
        else:
            hi = mid
    return hi


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

    _, _, _, _, p0, beta, zeta, segments, first = build(args)
    pvals = np.linalg.eigvalsh(sym(p0))
    mstar = threshold_m(p0, beta, zeta, args.tol)
    low_m, vals_m = min_eig_for_m(p0, beta, zeta, args.M)
    print(
        f"endpoint P0 finite-rank bound lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    print(
        f"  P0 min={pvals[0]:.12e} second={pvals[1]:.12e} "
        f"neg={(pvals < -args.tol).sum():3d}"
    )
    print(f"  scalar threshold M ~= {mstar:.12e}")
    print(
        f"  with M={args.M:g}: min={low_m:.12e} "
        f"neg={(vals_m < -args.tol).sum():3d}"
    )
    print("  low corrected eigenvalues: " + " ".join(f"{v:.3e}" for v in vals_m[:8]))


if __name__ == "__main__":
    main()
