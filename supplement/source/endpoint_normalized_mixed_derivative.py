#!/usr/bin/env python3
"""Mixed-derivative tests for the endpoint-normalized kernel.

For

  Khat(s,t)=P0(c e^s,c e^t)/[(c e^s-c)(c e^t-c)],

test whether partial_s partial_t Khat is positive definite.  Positivity would
give a direct double-Volterra proof of conditional positivity after endpoint
terms; failure means the normalized theorem still needs the coupled row
correction and not just a mixed-derivative kernel.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_normalized_mixed_derivative.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def khat(c: float, s: float, t: float) -> float:
    lam = c * math.exp(s)
    mu = c * math.exp(t)
    return p0_kernel(c, lam, mu) / ((lam - c) * (mu - c))


def mixed(c: float, s: float, t: float, h: float) -> float:
    return (
        khat(c, s + h, t + h)
        - khat(c, s + h, t - h)
        - khat(c, s - h, t + h)
        + khat(c, s - h, t - h)
    ) / (4.0 * h * h)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=70)
    parser.add_argument("--h", type=float, default=1e-4)
    parser.add_argument("--start", type=float, default=1e-3)
    parser.add_argument("--tol", type=float, default=1e-8)
    args = parser.parse_args()

    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = args.start + 0.5 * (args.T - args.start) * (nodes + 1.0)
    s_wts = 0.5 * (args.T - args.start) * weights
    root = np.sqrt(s_wts)
    mat = np.zeros((args.order, args.order), dtype=float)
    for i, s in enumerate(s_pts):
        for j, t in enumerate(s_pts[: i + 1]):
            val = root[i] * mixed(args.c, float(s), float(t), args.h) * root[j]
            mat[i, j] = mat[j, i] = val
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"endpoint normalized mixed derivative s=[{args.start:g},{args.T:g}] "
        f"order={args.order} h={args.h:g}"
    )
    print(
        f"  min={vals[0]:.12e} second={vals[1]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -args.tol).sum():3d}"
    )
    print("  low eigenvalues: " + " ".join(f"{x:.3e}" for x in vals[:10]))


if __name__ == "__main__":
    main()
