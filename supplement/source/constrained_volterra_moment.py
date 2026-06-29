#!/usr/bin/env python3
"""Galerkin test of the Volterra moment form on the endpoint-defect complement.

Represent f on [0,L] in the orthonormal shifted Legendre basis and impose
sampled endpoint-jet constraints

  Lambda_a(f) = sum_{k=0}^8 e_k(a) f^(k)(a)/k! = 0

for active endpoint centers a < s_*.  The script then projects the full
finite-core K_red Galerkin matrix to the numerical nullspace of those
constraints and reports positivity there.

This is a finite-dimensional diagnostic for the proposed smooth complement

  M_smooth = { f : Lambda_a(f)=0 for active endpoint defects a }.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402
from numpy.polynomial.legendre import Legendre  # noqa: E402

from endpoint_defect_family_mp import defect_at  # noqa: E402
from legendre_certificate import project_to_legendre, weighted_kernel  # noqa: E402


def shifted_legendre_derivative(order, length, x, deriv):
    z = 2.0 * x / length - 1.0
    out = np.zeros(order, dtype=float)
    scale = (2.0 / length) ** deriv
    for n in range(order):
        poly = Legendre.basis(n).deriv(deriv)
        out[n] = math.sqrt((2 * n + 1) / length) * scale * poly(z)
    return out


def lambda_row(center, e_coeffs, basis_order, length):
    row = np.zeros(basis_order, dtype=float)
    for k, coeff in enumerate(e_coeffs):
        deriv = shifted_legendre_derivative(basis_order, length, center, k)
        row += float(coeff) * deriv / math.factorial(k)
    return row


def sampled_centers(lo, hi, samples):
    if samples == 1:
        return [0.5 * (lo + hi)]
    return list(np.linspace(lo, hi, samples))


def constraint_matrix(args):
    rows = []
    centers = sampled_centers(args.constraint_min, args.constraint_max, args.constraints)
    for center in centers:
        vals, neg, vec = defect_at(
            mp.mpf(str(center)),
            args.jet_order,
            mp.mpf(str(args.endpoint_rmax)),
            args.endpoint_order,
            mp.mpf(str(args.endpoint_tol)),
        )
        if vals[0] >= 0:
            print(
                f"  warning: constraint center {center:.8g} has nonnegative "
                f"lambda0={float(vals[0]):.12e}"
            )
        rows.append(lambda_row(center, vec, args.basis, args.L))
    return centers, np.vstack(rows)


def nullspace(rows, tol):
    if rows.size == 0:
        return np.eye(rows.shape[1]), np.array([])
    u, s, vh = np.linalg.svd(rows, full_matrices=True)
    rank = int((s > tol * max(1.0, s[0] if len(s) else 1.0)).sum())
    return vh[rank:].T, s


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--basis", type=int, default=28)
    parser.add_argument("--quad", type=int, default=140)
    parser.add_argument("--laguerre", type=int, default=140)
    parser.add_argument("--constraints", type=int, default=8)
    parser.add_argument("--constraint-min", type=float, default=0.05)
    parser.add_argument("--constraint-max", type=float, default=0.54)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=50)
    parser.add_argument("--endpoint-rmax", type=float, default=12.0)
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--svd-tol", type=float, default=1e-10)
    args = parser.parse_args()

    pts, wts, _, full_w = weighted_kernel(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    kmat = project_to_legendre(full_w, pts, wts, args.basis, args.L)
    kmat = 0.5 * (kmat + kmat.T)

    centers, rows = constraint_matrix(args)
    ns, singular = nullspace(rows, args.svd_tol)
    constrained = ns.T @ kmat @ ns
    unconstrained_vals = np.linalg.eigvalsh(kmat)
    constrained_vals = (
        np.linalg.eigvalsh(0.5 * (constrained + constrained.T))
        if ns.shape[1]
        else np.array([])
    )
    row_resid = np.linalg.norm(rows @ ns) if ns.shape[1] else 0.0

    print(
        f"constrained Volterra moment kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} basis={args.basis} quad={args.quad} laguerre={args.laguerre}"
    )
    print(
        f"  constraints={args.constraints} interval=[{args.constraint_min:g},"
        f"{args.constraint_max:g}] jet_order={args.jet_order}"
    )
    print("  centers=" + ", ".join(f"{c:.6g}" for c in centers))
    print(
        f"  constraint rank={args.basis - ns.shape[1]} nullity={ns.shape[1]} "
        f"row_resid={row_resid:.12e}"
    )
    if len(singular):
        print(
            "  singular values="
            + " ".join(f"{value:.3e}" for value in singular[: min(12, len(singular))])
        )
    print(
        f"  unconstrained K min={unconstrained_vals[0]:.12e} "
        f"second={unconstrained_vals[1]:.12e}"
    )
    if len(constrained_vals):
        print(
            f"  constrained K min={constrained_vals[0]:.12e} "
            f"second={constrained_vals[1] if len(constrained_vals) > 1 else float('nan'):.12e}"
        )
        print(
            "  constrained low="
            + " ".join(f"{value:.3e}" for value in constrained_vals[: min(8, len(constrained_vals))])
        )
    else:
        print("  constrained space is zero-dimensional")


if __name__ == "__main__":
    main()
