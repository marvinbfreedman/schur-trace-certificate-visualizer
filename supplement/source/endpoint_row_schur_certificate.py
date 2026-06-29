#!/usr/bin/env python3
"""Two-row Schur certificate for the P0 and C Hardy completions.

Let A be P0 or C, and let R have rows beta,zeta in weighted coordinates.
The quotient theorem is A>=0 on ker R.  In finite dimension this is equivalent
to a singular Schur statement:

  A22 >= 0,  range(A12) subset range(A22),
  S = A11 - A12 A22^+ A12^T

relative to an orthonormal split row(R)^* plus ker(R).  A correction

  A + R^T M R

is positive iff the row-coordinate correction dominates -S.  This script
computes the intrinsic 2x2 Schur matrix and the minimal row-space correction.
It is a sharper target than the scalar bound M=14.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_certificate import build, sym

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_row_schur_certificate.py requires numpy") from exc


def split_schur(mat, rows, tol):
    # rows is k x n.  V spans the row space of rows, W its kernel.
    u, singular, vt = np.linalg.svd(rows, full_matrices=True)
    k = int((singular > 1e-12).sum())
    v = vt[:k].T
    w = vt[k:].T
    mat = sym(mat)
    a11 = v.T @ mat @ v
    a12 = v.T @ mat @ w
    a22 = sym(w.T @ mat @ w)
    vals22, vecs22 = np.linalg.eigh(a22)
    keep = vals22 > tol
    inv22 = (vecs22[:, keep] / vals22[keep]) @ vecs22[:, keep].T
    schur = sym(a11 - a12 @ inv22 @ a12.T)
    defect = 0.0
    if (~keep).any():
        defect = float(np.linalg.norm(a12 @ vecs22[:, ~keep]))

    # Convert the row-coordinate correction Gamma=-schur into the original
    # beta,zeta row metric: R^T M R = V Gamma V^T.  This matrix is diagnostic
    # only; the scalar bound M I is tested by S + M Sigma^2 >= 0 below.
    sigma = np.diag(singular[:k])
    sigma_inv = np.diag(1.0 / singular[:k])
    needed_gamma = sym(-schur)
    needed_m = sym(u[:, :k] @ sigma_inv @ needed_gamma @ sigma_inv @ u[:, :k].T)
    scaled = sym(sigma_inv @ (-schur) @ sigma_inv)
    scalar_threshold = max(0.0, float(np.linalg.eigvalsh(scaled)[-1]))
    return {
        "singular": singular,
        "quotient_vals": vals22,
        "schur": schur,
        "defect": defect,
        "needed_m": needed_m,
        "scalar_threshold": scalar_threshold,
        "scalar_schur": sym(schur + scalar_threshold * sigma @ sigma),
    }


def eig_summary(vals, tol):
    return f"min={vals[0]: .12e} neg={(vals < -tol).sum():3d}"


def report(name, mat, rows, scalar_m, tol):
    out = split_schur(mat, rows, tol)
    vals = np.linalg.eigvalsh(sym(mat))
    qvals = out["quotient_vals"]
    svals = np.linalg.eigvalsh(out["schur"])
    mvals = np.linalg.eigvalsh(sym(out["needed_m"]))
    scalar_svals = np.linalg.eigvalsh(out["scalar_schur"])
    corrected = sym(mat + scalar_m * rows.T @ rows)
    cvals = np.linalg.eigvalsh(corrected)
    print(f"  {name}:")
    print(f"    full      {eig_summary(vals, tol)}")
    print(f"    quotient  {eig_summary(qvals, tol)}")
    print(f"    Schur eig {' '.join(f'{x:.12e}' for x in svals)}")
    print(f"    range defect = {out['defect']:.12e}")
    print("    row singular values: " + " ".join(f"{x:.12e}" for x in out["singular"]))
    print(f"    Schur scalar threshold M ~= {out['scalar_threshold']:.12e}")
    print(
        "    threshold Schur eig "
        + " ".join(f"{x:.12e}" for x in scalar_svals)
    )
    print("    needed M-row matrix:")
    for row in out["needed_m"]:
        print("      " + " ".join(f"{x: .12e}" for x in row))
    print(f"    needed M eig {' '.join(f'{x:.12e}' for x in mvals)}")
    print(f"    scalar {scalar_m:g} correction {eig_summary(cvals, tol)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--M", type=float, default=14.0)
    args = parser.parse_args()

    p0, cmat, beta, zeta = build(args)
    rows = np.vstack([beta, zeta])
    print(
        f"endpoint row Schur certificate lambda=[c,c exp({args.T:g})] "
        f"order={args.order}"
    )
    report("P0", p0, rows, args.M, args.tol)
    report("C", cmat, rows, args.M, args.tol)


if __name__ == "__main__":
    main()
