#!/usr/bin/env python3
"""Endpoint-zero normalization of the beta/zeta quotient rows.

Since u_lambda vanishes at lambda=c, zeta(lambda)=1-c/lambda is a natural
positive endpoint factor.  Congruence by diag(1/zeta) preserves inertia away
from the endpoint and turns the two quotient rows into

  1,  b(lambda)=beta(lambda)/zeta(lambda).

This script samples b and the normalized P0/C kernels, looking for monotonicity
or a simpler row structure that could support a two-negative-square proof.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import (
    beta_closed,
    c_kernel,
    p0_kernel,
    zeta_closed,
)

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_row_normalization.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def finite_diff(values, xs):
    return np.diff(values) / np.diff(xs)


def build_kernel(args, kernel):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s = 0.5 * args.T * (nodes + 1.0)
    w = 0.5 * args.T * weights
    root = np.sqrt(w)
    lam = args.c * np.exp(s)
    z = np.array([zeta_closed(args.c, float(x)) for x in lam])
    mat = np.zeros((args.order, args.order), dtype=float)
    for i, x in enumerate(lam):
        for j, y in enumerate(lam[: i + 1]):
            val = root[i] * kernel(args.c, float(x), float(y)) * root[j] / (
                z[i] * z[j]
            )
            mat[i, j] = mat[j, i] = val
    b = root * np.array([beta_closed(args.c, float(x)) / zeta_closed(args.c, float(x)) for x in lam])
    one = root.copy()
    return s, lam, sym(mat), one, b


def build_lambda_minus_c_kernel(args, kernel):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s = 0.5 * args.T * (nodes + 1.0)
    w = 0.5 * args.T * weights
    root = np.sqrt(w)
    lam = args.c * np.exp(s)
    endpoint = lam - args.c
    mat = np.zeros((args.order, args.order), dtype=float)
    for i, x in enumerate(lam):
        for j, y in enumerate(lam[: i + 1]):
            val = root[i] * kernel(args.c, float(x), float(y)) * root[j] / (
                endpoint[i] * endpoint[j]
            )
            mat[i, j] = mat[j, i] = val
    exp_row = root / lam
    d_row = root * np.array(
        [beta_closed(args.c, float(x)) / (float(x) - args.c) for x in lam]
    )
    return s, lam, sym(mat), exp_row, d_row


def eig_line(name, mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<16} min={vals[0]: .12e} second={vals[1]: .12e} "
        f"neg={(vals < -tol).sum():3d}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=100)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    s_grid = np.linspace(1e-6, args.T, args.samples)
    lam_grid = args.c * np.exp(s_grid)
    beta = np.array([beta_closed(args.c, float(x)) for x in lam_grid])
    zeta = np.array([zeta_closed(args.c, float(x)) for x in lam_grid])
    b = beta / zeta
    db = finite_diff(b, s_grid)
    d2b = finite_diff(db, 0.5 * (s_grid[1:] + s_grid[:-1]))

    print(f"endpoint row normalization samples s in [0,{args.T:g}]")
    print(f"  b(0+) approx={b[0]:.12e} b(T)={b[-1]:.12e}")
    print(f"  db min={db.min():.12e} max={db.max():.12e}")
    print(f"  d2b min={d2b.min():.12e} max={d2b.max():.12e}")
    print("  sample b:")
    for idx in np.linspace(0, args.samples - 1, min(8, args.samples), dtype=int):
        print(f"    s={s_grid[idx]:.6g} b={b[idx]:.12e}")

    for name, kernel in (("P0/zeta", p0_kernel), ("C/zeta", c_kernel)):
        _, _, mat, one, brow = build_kernel(args, kernel)
        eig_line(name, mat, args.tol)
        corrected = mat + 14.0 * (np.outer(one, one) + np.outer(brow, brow))
        eig_line(name + "+rows", corrected, args.tol)

    for name, kernel in (("P0/(lam-c)", p0_kernel), ("C/(lam-c)", c_kernel)):
        _, _, mat, erow, drow = build_lambda_minus_c_kernel(args, kernel)
        eig_line(name, mat, args.tol)
        corrected = mat + 14.0 * (np.outer(erow, erow) + np.outer(drow, drow))
        eig_line(name + "+rows", corrected, args.tol)


if __name__ == "__main__":
    main()
