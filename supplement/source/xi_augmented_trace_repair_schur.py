#!/usr/bin/env python3
r"""Finite Schur certificate for the augmented trace repair (Lambda, Mu).

The Mellin-boundary concomitant adds rows

    Mu_z(f)=int_0^L B(s,z)f(s)ds

to the old endpoint trace matrix R_Lambda.  This script asks whether the
augmented quotient repair can be chosen positive:

    P_aug = K + R_aug^* D_aug R_aug >= 0,     D_aug >= 0,

and whether the boundary Mu rows are annihilated on ker R_aug.  This is the
finite model for the next analytic theorem: adding the Mellin-boundary trace
must not introduce a negative Schur defect.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_trace_map_constructor import ctranspose, to_float  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    diag,
    frob_norm,
    gram_matrix,
    max_eig_or_zero,
    minmax_eigs,
    positive_part_inverse,
    shifted_legendre_polys,
    trace_matrix,
)
from xi_mellin_boundary_concomitant import (  # noqa: E402
    local_trace_args,
    prefix_rows,
    split_trace_nullspace_complex,
    stack_rows,
)


def herm(a: mp.matrix) -> mp.matrix:
    return (a + ctranspose(a)) / 2


def complex_minmax_eigs(mat: mp.matrix):
    vals = mp.eighe(herm(mat), eigvals_only=True)
    return vals, vals[0], vals[-1]


def complex_op_norm(mat: mp.matrix):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eighe(herm(ctranspose(mat) * mat), eigvals_only=True)
    return mp.sqrt(max(mp.mpf("0"), vals[-1]))


def complex_columns(mat: mp.matrix, indices):
    return columns(mat, indices)


def split_for_repair(K: mp.matrix, R: mp.matrix, rank_tol_text: str):
    N, U, rvals, rank_tol = split_trace_nullspace_complex(R, rank_tol_text)
    return N, U, rvals, rank_tol


def positive_part_inverse_herm(mat: mp.matrix, tol: mp.mpf):
    vals, vecs = mp.eighe(herm(mat), eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > tol]
    zero = [i for i, val in enumerate(vals) if val <= tol]
    inv = mp.matrix(mat.rows)
    if keep:
        vkeep = complex_columns(vecs, keep)
        inv = vkeep * diag([1 / vals[i] for i in keep]) * ctranspose(vkeep)
    return vals, vecs, keep, zero, inv


def quotient_repair_certificate(K: mp.matrix, R: mp.matrix, args):
    N, U, rvals, rank_tol = split_for_repair(K, R, args.rank_tol)
    A = ctranspose(N) * K * N
    B = ctranspose(N) * K * U
    C = ctranspose(U) * K * U

    if A.rows:
        avals, avecs, _a_keep, a_zero, Aplus = positive_part_inverse_herm(A, mp.mpf(args.psd_tol))
    else:
        avals, avecs, a_zero, Aplus = [], mp.matrix(0), [], mp.matrix(0)
    range_resid = mp.mpf("0")
    if a_zero and B.cols:
        zvecs = complex_columns(avecs, a_zero)
        range_resid = frob_norm(ctranspose(zvecs) * B)
    gamma2 = ctranspose(B) * Aplus * B if B.rows and B.cols else mp.matrix(U.cols)
    H = gamma2 - C if B.cols else -C
    hvals, hmin, hmax = complex_minmax_eigs(H) if H.rows else ([], mp.mpf("0"), mp.mpf("0"))
    alpha = max(mp.mpf("0"), hmax) + mp.mpf(args.margin)
    M = alpha * mp.eye(U.cols)

    D = mp.matrix(R.rows)
    if U.cols:
        RU = R * U
        gram_ru = herm(ctranspose(RU) * RU)
        sig_vals, sig_vecs = mp.eighe(gram_ru, eigvals_only=False)
        keep = [i for i, val in enumerate(sig_vals) if val > mp.mpf(args.rank_tol)]
        if keep:
            V = complex_columns(sig_vecs, keep)
            inv_sqrt = V * diag([1 / mp.sqrt(sig_vals[i]) for i in keep]) * ctranspose(V)
            W = RU * inv_sqrt
            D = W * (alpha * mp.eye(W.cols)) * ctranspose(W)

    P_from_trace = K + ctranspose(R) * D * R
    pvals, pmin, pmax = complex_minmax_eigs(P_from_trace)
    dvals, dmin, dmax = complex_minmax_eigs(D) if D.rows else ([], mp.mpf("0"), mp.mpf("0"))
    kvals, kmin, kmax = complex_minmax_eigs(K)
    avals_full, amin, amax = complex_minmax_eigs(A) if A.rows else ([], mp.mpf("0"), mp.mpf("0"))
    b_norm = complex_op_norm(B)
    normalized_range = range_resid / b_norm if b_norm else mp.mpf("0")
    schur = C + M - gamma2 if U.cols else mp.matrix(0)
    svals, smin, smax = complex_minmax_eigs(schur) if U.cols else ([], mp.mpf("0"), mp.mpf("0"))
    return {
        "rank": U.cols,
        "nullity": N.cols,
        "rankTolerance": to_float(rank_tol),
        "kMin": to_float(kmin),
        "kMax": to_float(kmax),
        "aMin": to_float(amin),
        "aMax": to_float(amax),
        "crossNorm": to_float(b_norm),
        "rangeResidual": to_float(range_resid),
        "normalizedRangeResidual": to_float(normalized_range),
        "schurDefectMin": to_float(hmin),
        "schurDefectMax": to_float(hmax),
        "repairAlpha": to_float(alpha),
        "constructedSchurMin": to_float(smin),
        "constructedSchurMax": to_float(smax),
        "dMin": to_float(dmin),
        "dMax": to_float(dmax),
        "pMin": to_float(pmin),
        "pMax": to_float(pmax),
        "lowK": [to_float(v) for v in kvals[: min(6, len(kvals))]],
        "lowA": [to_float(v) for v in avals_full[: min(6, len(avals_full))]],
        "lowSchurDefect": [to_float(v) for v in hvals[: min(6, len(hvals))]],
        "lowP": [to_float(v) for v in pvals[: min(6, len(pvals))]],
        "lowD": [to_float(v) for v in dvals[: min(6, len(dvals))]] if D.rows else [],
        "positiveRepair": bool(dmin >= -mp.mpf(args.eig_tol) and pmin >= -mp.mpf(args.eig_tol)),
    }


def row_action_on_null(R_mu: mp.matrix, N: mp.matrix):
    if N.cols == 0:
        return mp.mpf("0")
    return complex_op_norm(R_mu * N)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="full")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=10)
    parser.add_argument("--quad", type=int, default=18)
    parser.add_argument("--constraints", type=int, default=7)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=20)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-24")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--eig-tol", default="1e-20")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--fourier-scale", default="0.5")
    parser.add_argument("--amp-kind", choices=["positive", "even"], default="positive")
    parser.add_argument("--hardy-normalization", action="store_true", default=True)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="xi_augmented_trace_repair_schur.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    L = mp.mpf(args.L)
    K, _ = gram_matrix(args, mp.mpf(args.omega), L)
    polys = shifted_legendre_polys(args.basis, L)
    _centers, R_lambda = trace_matrix(polys, local_trace_args(args))
    R_mu = prefix_rows(args, polys)
    R_aug = stack_rows(R_lambda, R_mu)

    old = quotient_repair_certificate(K, R_lambda, args)
    augmented = quotient_repair_certificate(K, R_aug, args)
    N_old, _U_old, _rv_old, _tol_old = split_for_repair(K, R_lambda, args.rank_tol)
    N_aug, _U_aug, _rv_aug, _tol_aug = split_for_repair(K, R_aug, args.rank_tol)
    old_mu_action = row_action_on_null(R_mu, N_old)
    aug_mu_action = row_action_on_null(R_mu, N_aug)
    data = {
        "theoremName": "augmented Mellin-boundary trace repair Schur certificate",
        "kind": args.kind,
        "omega": float(args.omega),
        "L": float(args.L),
        "basis": args.basis,
        "constraints": args.constraints,
        "lambdaRows": R_lambda.rows,
        "muRows": R_mu.rows,
        "augmentedRows": R_aug.rows,
        "oldLambdaRepair": old,
        "augmentedRepair": augmented,
        "oldMuActionOnLambdaNullspace": to_float(old_mu_action),
        "augmentedMuActionOnAugmentedNullspace": to_float(aug_mu_action),
        "augmentedRepairPositive": bool(augmented["positiveRepair"]),
        "augmentedMuAnnihilated": bool(aug_mu_action < mp.mpf(args.eig_tol)),
        "diagnosis": (
            "Finite Schur repair with R_aug=(Lambda,Mu).  If P=K+R_aug^*D R_aug "
            "is positive and D>=0 while Mu vanishes on ker R_aug, the augmented "
            "repair is finite-positive; the remaining task is the continuum "
            "closure/positivity proof."
        ),
        "nextProofTarget": (
            "Lift the finite augmented repair to the continuum: prove the Mu_z "
            "primitive rows are closed trace functionals, D_aug>=0 is bounded in "
            "the transported trace norm, and K+R_aug^*D_aug R_aug remains positive "
            "under Galerkin exhaustion."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented Mellin-boundary trace repair Schur certificate")
    print(
        f"  kind={args.kind} omega={args.omega} basis={args.basis} "
        f"constraints={args.constraints}"
    )
    print(f"  old repair pmin={old['pMin']:.3e} dmin={old['dMin']:.3e}")
    print(
        f"  augmented repair pmin={augmented['pMin']:.3e} "
        f"dmin={augmented['dMin']:.3e}"
    )
    print(f"  old Mu action on ker Lambda: {to_float(old_mu_action):.3e}")
    print(f"  augmented Mu action on ker augmented: {to_float(aug_mu_action):.3e}")
    print(f"  augmented repair positive: {data['augmentedRepairPositive']}")
    print(f"  augmented Mu annihilated: {data['augmentedMuAnnihilated']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
