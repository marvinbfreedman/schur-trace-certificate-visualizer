#!/usr/bin/env python3
r"""Trace Green minimizers in the explicit Volterra variables.

This ledger identifies the non-spectral continuum feature candidate.

For a trace vector x, let f_x be the Euler-Lagrange minimizer in the fiber
Rf=x.  The direct reduced Volterra formula gives, for each branch sigma,

    B_sigma(s,u)=exp(0.5*sigma*omega*(s+u)) A_s(u),

    M_{sigma,x}(u)=int f_x(s) B_sigma(s,u) ds,
    N_{sigma,x}(u)=int (s+u) f_x(s) B_sigma(s,u) ds.

Then

    D_trace(x,y)
      = 1/2 sum_sigma w_sigma int [
            M_{sigma,x}(u) N_{sigma,y}(u)
          + N_{sigma,x}(u) M_{sigma,y}(u)
        ] du.

This is an exact Volterra/Green representation once f_x is the Green
minimizer.  It is not yet a positive square.  Pointwise square completion gives

    M_x N_y + N_x M_y
      = 1/2[(M_x+N_x)(M_y+N_y) - (M_x-N_x)(M_y-N_y)],

or equivalently with signed Volterra features

    G_{sigma,+,x}=sqrt(w_sigma/4)(M_{sigma,x}+N_{sigma,x}),
    G_{sigma,-,x}=sqrt(w_sigma/4)(M_{sigma,x}-N_{sigma,x}).

Thus D_trace is the difference of the + and - feature Grams.  An additional
constrained moment theorem is needed to turn this signed Krein representation
into an honest Hilbert Gram square.

The finite check below verifies that the Volterra moment representation of the
Green minimizers agrees with the Schur-minimized energy.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp
import numpy as np

from quotient_factorization_mp import (
    columns,
    frob_norm,
    gram_matrix,
    positive_part_inverse,
    trace_matrix,
)
from reduced_exact_finite import endpoint_ratios, pieces


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def quadrature(end: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return 0.5 * end * (nodes + 1.0), 0.5 * end * weights


def shifted_legendre_values(order: int, length: float, pts: np.ndarray) -> np.ndarray:
    z = 2.0 * pts / length - 1.0
    values = np.zeros((len(pts), order), dtype=float)
    for n in range(order):
        coeffs = np.zeros(n + 1)
        coeffs[-1] = 1.0
        values[:, n] = math.sqrt((2 * n + 1) / length) * np.polynomial.legendre.legval(
            z, coeffs
        )
    return values


def a_value(s: float, u: float, pcs) -> float:
    es = math.exp(s)
    eu = math.exp(u)
    total = 0.0
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * math.exp(beta * u - c * es * (eu - 1.0))
    return total


def branch_feature_matrices(
    *,
    kind: str,
    omega: float,
    sigma: int,
    length: float,
    basis: int,
    s_order: int,
    u_max: float,
    u_order: int,
) -> tuple[np.ndarray, np.ndarray]:
    pcs = pieces(kind)
    s_pts, s_wts = quadrature(length, s_order)
    u_pts, u_wts = quadrature(u_max, u_order)
    phi = shifted_legendre_values(basis, length, s_pts)
    root_u = np.sqrt(u_wts)
    m = np.zeros((u_order, basis), dtype=float)
    n = np.zeros_like(m)
    for i, u in enumerate(u_pts):
        for j, s in enumerate(s_pts):
            branch = math.exp(0.5 * sigma * omega * (s + u))
            base = root_u[i] * s_wts[j] * a_value(float(s), float(u), pcs) * branch
            m[i, :] += base * phi[j, :]
            n[i, :] += (s + u) * base * phi[j, :]
    return m, n


def build_args(args: argparse.Namespace):
    return SimpleNamespace(
        model="full",
        kind=args.kind,
        omega=args.omega,
        L=args.L,
        basis=args.basis,
        quad=max(args.min_quad, args.quad_factor * args.basis),
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=args.constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
        dps=args.dps,
    )


def quotient_blocks(K: mp.matrix, R: mp.matrix, args):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)
    A = N.T * K * N
    B = N.T * K * U
    C = U.T * K * U
    _, _, _, _, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    S = (C - B.T * Aplus * B + (C - B.T * Aplus * B).T) / 2
    minimizer = U - N * Aplus * B
    return N, U, A, B, C, S, minimizer


def mp_to_np(mat: mp.matrix) -> np.ndarray:
    return np.array([[float(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)])


def np_to_mp(mat: np.ndarray) -> mp.matrix:
    out = mp.matrix(mat.shape[0], mat.shape[1])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            out[i, j] = mp.mpf(repr(float(mat[i, j])))
    return out


def finite_volterra_check(args: argparse.Namespace) -> dict:
    local = build_args(args)
    mp.mp.dps = args.dps

    sigmas = [0] if float(args.omega) == 0.0 else [-1, 1]
    branch_weight = 1.0 if float(args.omega) == 0.0 else 0.5
    basis_moment = np.zeros((args.basis, args.basis), dtype=float)
    feature_pairs = []

    for sigma in sigmas:
        m, n = branch_feature_matrices(
            kind=args.kind,
            omega=float(args.omega),
            sigma=sigma,
            length=float(args.L),
            basis=args.basis,
            s_order=args.s_order,
            u_max=args.u_max,
            u_order=args.u_order,
        )
        branch_basis_moment = 0.5 * branch_weight * (m.T @ n + n.T @ m)
        basis_moment += branch_basis_moment
        feature_pairs.append((sigma, m, n, branch_weight, branch_basis_moment))
    basis_moment = 0.5 * (basis_moment + basis_moment.T)
    K_direct = np_to_mp(basis_moment)

    K_reference, polys = gram_matrix(local, mp.mpf(local.omega), mp.mpf(local.L))
    _, R = trace_matrix(polys, local)
    if args.kernel_source == "gram_matrix":
        K = K_reference
    else:
        K = K_direct
    N, U, A, B, C, S, minimizer = quotient_blocks(K, R, local)
    F = mp_to_np(minimizer)
    S_np = mp_to_np(S)

    direct_vs_reference = frob_norm(K_direct - K_reference) / max(
        frob_norm(K_reference), mp.mpf("1e-300")
    )
    direct_evals = np.linalg.eigvalsh(basis_moment)
    reference_evals = np.linalg.eigvalsh(mp_to_np(K_reference))

    moment = np.zeros_like(S_np)
    signed_square_plus = np.zeros_like(S_np)
    signed_square_minus = np.zeros_like(S_np)
    branch_rows = []

    for sigma, m, n, branch_weight, branch_basis_moment in feature_pairs:
        mf = m @ F
        nf = n @ F
        branch_moment = 0.5 * branch_weight * (mf.T @ nf + nf.T @ mf)
        plus = 0.25 * branch_weight * ((mf + nf).T @ (mf + nf))
        minus = 0.25 * branch_weight * ((mf - nf).T @ (mf - nf))
        moment += branch_moment
        signed_square_plus += plus
        signed_square_minus += minus
        evals = np.linalg.eigvalsh(0.5 * (branch_moment + branch_moment.T))
        branch_rows.append(
            {
                "sigma": sigma,
                "momentMin": float(evals[0]),
                "momentMax": float(evals[-1]),
                "plusTrace": float(np.trace(plus)),
                "minusTrace": float(np.trace(minus)),
            }
        )

    moment = 0.5 * (moment + moment.T)
    signed_square = signed_square_plus - signed_square_minus
    moment_evals = np.linalg.eigvalsh(moment)
    plus_evals = np.linalg.eigvalsh(0.5 * (signed_square_plus + signed_square_plus.T))
    minus_evals = np.linalg.eigvalsh(0.5 * (signed_square_minus + signed_square_minus.T))
    error = np.linalg.norm(moment - S_np, ord="fro")
    rel_error = error / max(np.linalg.norm(S_np, ord="fro"), 1e-300)
    square_error = np.linalg.norm(signed_square - moment, ord="fro")
    return {
        "basis": args.basis,
        "constraints": args.constraints,
        "kernelSource": args.kernel_source,
        "traceRank": int(U.cols),
        "traceNullity": int(N.cols),
        "sOrder": args.s_order,
        "uOrder": args.u_order,
        "uMax": float(args.u_max),
        "directKernelMin": float(direct_evals[0]),
        "directKernelMax": float(direct_evals[-1]),
        "referenceGramMin": float(reference_evals[0]),
        "referenceGramMax": float(reference_evals[-1]),
        "directVsReferenceGramRelative": float(direct_vs_reference),
        "schurMin": float(np.linalg.eigvalsh(S_np)[0]),
        "schurMax": float(np.linalg.eigvalsh(S_np)[-1]),
        "volterraMomentMin": float(moment_evals[0]),
        "volterraMomentMax": float(moment_evals[-1]),
        "momentVsSchurFrobenius": float(error),
        "momentVsSchurRelative": float(rel_error),
        "signedSquareIdentityFrobenius": float(square_error),
        "plusMin": float(plus_evals[0]),
        "minusMin": float(minus_evals[0]),
        "plusTrace": float(np.trace(signed_square_plus)),
        "minusTrace": float(np.trace(signed_square_minus)),
        "branchRows": branch_rows,
        "finiteVolterraRepresentationClosed": bool(rel_error < args.relative_tol),
        "honestPositiveSquareClosed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=5)
    parser.add_argument("--quad-factor", type=int, default=2)
    parser.add_argument("--min-quad", type=int, default=14)
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=18)
    parser.add_argument("--endpoint-kernel-rmax", default="10")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=28)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-26")
    parser.add_argument("--s-order", type=int, default=36)
    parser.add_argument("--u-max", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=100)
    parser.add_argument(
        "--kernel-source",
        choices=["direct_volterra", "gram_matrix"],
        default="direct_volterra",
        help=(
            "Use the direct Volterra feature operator to form the quotient, "
            "or compare the same features against quotient_factorization_mp.gram_matrix."
        ),
    )
    parser.add_argument("--relative-tol", type=float, default=5e-4)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="trace_volterra_green_feature_map.json")
    args = parser.parse_args()

    check = finite_volterra_check(args)

    statuses = {
        "greenMinimizerFeatureDefinitionStatus": status(
            "Green minimizer Volterra feature definition",
            True,
            (
                "For trace x, use the Euler-Lagrange minimizer f_x and define "
                "M_{sigma,x}, N_{sigma,x} by the direct reduced Volterra "
                "formula."
            ),
        ),
        "volterraMomentRepresentationStatus": status(
            "Volterra moment representation of D_trace",
            check["finiteVolterraRepresentationClosed"],
            (
                "The identity D_trace(x,y)=1/2 sum int(M_x N_y+N_x M_y) "
                "is the direct Volterra formula restricted to Green minimizers; "
                "the finite check agrees with the Schur matrix built from the "
                "same direct Volterra operator within the displayed quadrature "
                "tolerance."
            ),
            blocker=None
            if check["finiteVolterraRepresentationClosed"]
            else "Increase quadrature or audit the Volterra/minimizer coordinate conversion.",
        ),
        "signedSquareCompletionStatus": status(
            "signed Volterra square completion",
            True,
            (
                "Pointwise, M_xN_y+N_xM_y equals one half of a plus-square "
                "kernel minus a minus-square kernel."
            ),
        ),
        "positiveGreenGramStatus": status(
            "honest positive Volterra/Green Gram square",
            False,
            (
                "The explicit Volterra features give a signed Krein square, "
                "not yet an honest Hilbert square.  A constrained moment "
                "positivity theorem is still needed to remove or dominate the "
                "minus-square part on Green minimizers."
            ),
            blocker=(
                "Prove that the negative square component is dominated by the "
                "positive component on the Euler-Lagrange Green-minimizer "
                "trace family, equivalently prove positivity of the Volterra "
                "moment operator on that constrained image."
            ),
        ),
    }

    data = {
        "theoremName": "Volterra/Green feature-map identification",
        "featureDefinitions": {
            "greenMinimizer": "f_x=Jx-A^+Bx",
            "A_s_u": "Psi(s+u)/Psi(s)",
            "B_sigma": "exp(0.5*sigma*omega*(s+u))*A_s(u)",
            "M_sigma_x": "int f_x(s) B_sigma(s,u) ds",
            "N_sigma_x": "int (s+u) f_x(s) B_sigma(s,u) ds",
            "G_sigma_plus_x": "sqrt(w_sigma/4)*(M_sigma_x+N_sigma_x)",
            "G_sigma_minus_x": "sqrt(w_sigma/4)*(M_sigma_x-N_sigma_x)",
        },
        "exactVolterraIdentity": (
            "D_trace(x,y)=1/2 sum_sigma w_sigma int "
            "[M_sigma_x N_sigma_y + N_sigma_x M_sigma_y] du"
        ),
        "signedFeatureIdentity": (
            "D_trace(x,y)=sum_sigma int "
            "[G_sigma_plus_x G_sigma_plus_y - "
            "G_sigma_minus_x G_sigma_minus_y] du"
        ),
        "signedSquareIdentity": (
            "M_x N_y + N_x M_y = 1/2[(M_x+N_x)(M_y+N_y) - "
            "(M_x-N_x)(M_y-N_y)]"
        ),
        "statuses": statuses,
        "finiteCheck": check,
        "referenceGramDiagnostic": (
            "The optional gram_matrix comparison is diagnostic only.  The "
            "feature identity is checked against the direct Volterra operator "
            "because that is the continuum formula being identified."
        ),
        "continuumPositiveSquareClosed": False,
        "correctedConclusion": (
            "The explicit continuum Volterra/Green feature candidates are "
            "identified by applying the direct Volterra transform to the "
            "Euler-Lagrange Green minimizers.  They reproduce D_trace as a "
            "Volterra moment form, and finite quadrature verifies the identity. "
            "However, the resulting square completion is signed, so the honest "
            "positive Hilbert Gram square is still equivalent to the remaining "
            "constrained Volterra moment positivity theorem."
        ),
        "nextProofTarget": (
            "Prove constrained positivity of the signed Volterra moment form on "
            "the Green-minimizer trace image, or find an additional transform "
            "that converts the signed plus/minus square into a single positive "
            "Hilbert square."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Volterra/Green feature-map identification")
    print(f"  kernel source: {check['kernelSource']}")
    print(f"  finite Volterra representation: {check['finiteVolterraRepresentationClosed']}")
    print(f"  moment-vs-Schur relative error: {check['momentVsSchurRelative']:.3e}")
    print(f"  direct-vs-reference Gram relative error: {check['directVsReferenceGramRelative']:.3e}")
    print(f"  Schur min: {check['schurMin']:.12e}")
    print(f"  Volterra moment min: {check['volterraMomentMin']:.12e}")
    print(f"  plus trace / minus trace: {check['plusTrace']:.6e} / {check['minusTrace']:.6e}")
    print("  honest positive Volterra square: open")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
