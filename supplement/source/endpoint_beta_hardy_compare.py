#!/usr/bin/env python3
"""Compare the endpoint beta-derivative form with the free Hardy square."""

from __future__ import annotations

import argparse

from endpoint_beta_derivative import build, sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_beta_hardy_compare.py requires numpy") from exc


def eig_line(name, mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<22} min={vals[0]: .12e} max={vals[-1]: .12e} "
        f"neg={(vals < -tol).sum():3d}"
    )
    return vals


def generalized_ratios(a, b, tol):
    vals, vecs = np.linalg.eigh(sym(b))
    keep = vals > tol
    if not keep.any():
        return None
    q = vecs[:, keep] / np.sqrt(vals[keep])
    return np.linalg.eigvalsh(sym(q.T @ sym(a) @ q))


def restrict_boundary_null(mat, boundary, tol):
    vals, vecs = np.linalg.eigh(sym(boundary))
    basis = vecs[:, vals <= tol]
    return sym(basis.T @ sym(mat) @ basis)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=np.pi)
    parser.add_argument("--beta", type=float, default=1.5)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--t-order", type=int, default=90)
    parser.add_argument("--rmax", type=float, default=24.0)
    parser.add_argument("--r-order", type=int, default=400)
    parser.add_argument("--tol", type=float, default=1e-10)
    args = parser.parse_args()

    mats = build(args)
    rform = mats["D_beta-.5N"]
    hardy = mats["free_Hardy"]
    boundary = mats["Boundary"]
    q2 = mats["2Q"]
    print(
        f"endpoint beta/Hardy comparison beta={args.beta:g} "
        f"T={args.T:g} order={args.t_order}"
    )
    eig_line("R=D_beta-.5N", rform, args.tol)
    eig_line("free Hardy", hardy, args.tol)
    eig_line("R-freeHardy", rform - hardy, args.tol)
    eig_line("2Q", q2, args.tol)

    rn = restrict_boundary_null(rform, boundary, args.tol)
    hn = restrict_boundary_null(hardy, boundary, args.tol)
    dn = restrict_boundary_null(rform - hardy, boundary, args.tol)
    eig_line("R on boundary-null", rn, args.tol)
    eig_line("Hardy on null", hn, args.tol)
    eig_line("R-Hardy on null", dn, args.tol)

    ratios = generalized_ratios(rn, hn, args.tol)
    if ratios is not None:
        print(
            "  R/freeHardy on boundary-null: "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e} dim={len(ratios)}"
        )


if __name__ == "__main__":
    main()
