#!/usr/bin/env python3
r"""Contraction test for the signed Volterra/Green feature split.

The trace-side feature identity gives

    D_trace = P - M,

where P is the plus-feature Gram and M is the minus-feature Gram on the
Euler-Lagrange Green-minimizer trace image.  This script computes the finite
generalized spectrum of

    M <= theta P

in the direct Volterra quotient model.  The target theorem is theta <= 1 in
the continuum limit.  A strict finite gap suggests a Hardy/Volterra
contraction; theta close to 1 suggests an exact boundary-null identity plus a
nonnegative residual square.
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

from quotient_factorization_mp import shifted_legendre_polys, trace_matrix
from trace_volterra_green_feature_map import (
    branch_feature_matrices,
    mp_to_np,
    np_to_mp,
    quotient_blocks,
    status,
)


def sym(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat + mat.T)


def generalized_psd_spectrum(numer: np.ndarray, denom: np.ndarray, tol: float):
    dvals, dvecs = np.linalg.eigh(sym(denom))
    cutoff = tol * max(1.0, float(np.max(np.abs(dvals))) if dvals.size else 0.0)
    keep = np.array([i for i, val in enumerate(dvals) if val > cutoff], dtype=int)
    null = np.array([i for i, val in enumerate(dvals) if val <= cutoff], dtype=int)
    if keep.size == 0:
        return {
            "values": np.array([], dtype=float),
            "vectors": np.zeros((denom.shape[0], 0), dtype=float),
            "denomEigenvalues": dvals,
            "keep": keep,
            "null": null,
            "cutoff": cutoff,
            "nullNumeratorMax": None,
        }
    vkeep = dvecs[:, keep]
    inv_root = vkeep @ np.diag(1.0 / np.sqrt(dvals[keep]))
    reduced = sym(inv_root.T @ sym(numer) @ inv_root)
    vals, vecs_red = np.linalg.eigh(reduced)
    vecs = inv_root @ vecs_red
    null_num = None
    if null.size:
        vnull = dvecs[:, null]
        null_block = sym(vnull.T @ sym(numer) @ vnull)
        null_num = float(np.linalg.eigvalsh(null_block)[-1])
    return {
        "values": vals,
        "vectors": vecs,
        "denomEigenvalues": dvals,
        "keep": keep,
        "null": null,
        "cutoff": cutoff,
        "nullNumeratorMax": null_num,
    }


def build_local_args(args: argparse.Namespace, basis: int, constraints: int):
    return SimpleNamespace(
        model="full",
        kind=args.kind,
        omega=args.omega,
        L=args.L,
        basis=basis,
        constraints=constraints,
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


def feature_grams(args: argparse.Namespace, basis: int):
    sigmas = [0] if float(args.omega) == 0.0 else [-1, 1]
    branch_weight = 1.0 if float(args.omega) == 0.0 else 0.5
    plus_basis = np.zeros((basis, basis), dtype=float)
    minus_basis = np.zeros_like(plus_basis)
    branch_rows = []
    for sigma in sigmas:
        m, n = branch_feature_matrices(
            kind=args.kind,
            omega=float(args.omega),
            sigma=sigma,
            length=float(args.L),
            basis=basis,
            s_order=args.s_order,
            u_max=args.u_max,
            u_order=args.u_order,
        )
        plus = 0.25 * branch_weight * ((m + n).T @ (m + n))
        minus = 0.25 * branch_weight * ((m - n).T @ (m - n))
        plus_basis += plus
        minus_basis += minus
        branch_rows.append(
            {
                "sigma": sigma,
                "plusTraceBasis": float(np.trace(plus)),
                "minusTraceBasis": float(np.trace(minus)),
                "plusMinBasis": float(np.linalg.eigvalsh(sym(plus))[0]),
                "minusMinBasis": float(np.linalg.eigvalsh(sym(minus))[0]),
            }
        )
    return sym(plus_basis), sym(minus_basis), branch_rows


def run_case(args: argparse.Namespace, basis: int, constraints: int) -> dict:
    local = build_local_args(args, basis, constraints)
    mp.mp.dps = args.dps
    plus_basis, minus_basis, branch_rows = feature_grams(args, basis)
    direct = sym(plus_basis - minus_basis)
    polys = shifted_legendre_polys(basis, mp.mpf(args.L))
    _, R = trace_matrix(polys, local)
    N, U, A, B, C, S, minimizer = quotient_blocks(np_to_mp(direct), R, local)
    F = mp_to_np(minimizer)
    plus = sym(F.T @ plus_basis @ F)
    minus = sym(F.T @ minus_basis @ F)
    defect = sym(plus - minus)
    schur = mp_to_np(S)
    spec = generalized_psd_spectrum(minus, plus, args.generalized_tol)
    values = spec["values"]
    top = float(values[-1]) if values.size else float("nan")
    gap = 1.0 - top if values.size else float("nan")
    top_vec = spec["vectors"][:, -1] if values.size else np.zeros(F.shape[1])
    if top_vec.size:
        pnorm = math.sqrt(max(float(top_vec.T @ plus @ top_vec), 0.0))
        if pnorm:
            top_vec = top_vec / pnorm
    top_plus = float(top_vec.T @ plus @ top_vec) if top_vec.size else float("nan")
    top_minus = float(top_vec.T @ minus @ top_vec) if top_vec.size else float("nan")
    top_defect = float(top_vec.T @ defect @ top_vec) if top_vec.size else float("nan")
    defect_eigs = np.linalg.eigvalsh(defect)
    schur_eigs = np.linalg.eigvalsh(sym(schur))
    schur_error = np.linalg.norm(defect - sym(schur), ord="fro")
    schur_rel_error = schur_error / max(np.linalg.norm(sym(schur), ord="fro"), 1e-300)
    return {
        "basis": basis,
        "constraints": constraints,
        "traceRank": int(U.cols),
        "traceNullity": int(N.cols),
        "sOrder": args.s_order,
        "uOrder": args.u_order,
        "uMax": float(args.u_max),
        "plusMin": float(np.linalg.eigvalsh(plus)[0]),
        "plusMax": float(np.linalg.eigvalsh(plus)[-1]),
        "minusMin": float(np.linalg.eigvalsh(minus)[0]),
        "minusMax": float(np.linalg.eigvalsh(minus)[-1]),
        "defectMin": float(defect_eigs[0]),
        "defectMax": float(defect_eigs[-1]),
        "schurMin": float(schur_eigs[0]),
        "schurMax": float(schur_eigs[-1]),
        "defectVsSchurRelative": float(schur_rel_error),
        "generalizedCutoff": float(spec["cutoff"]),
        "generalizedEigenvalues": [float(v) for v in values],
        "topContractionEigenvalue": top,
        "topContractionGapToOne": gap,
        "plusNullity": int(spec["null"].size),
        "plusNullMinusMax": spec["nullNumeratorMax"],
        "topVectorInTraceCoordinates": [float(v) for v in top_vec],
        "topVectorPlusEnergy": top_plus,
        "topVectorMinusEnergy": top_minus,
        "topVectorDefectEnergy": top_defect,
        "branchRows": branch_rows,
        "finiteDominationClosed": bool(
            values.size > 0
            and top <= 1.0 + args.close_tol
            and (
                spec["nullNumeratorMax"] is None
                or spec["nullNumeratorMax"] <= args.close_tol
            )
        ),
    }


def parse_ints(text: str) -> list[int]:
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis-list", default="6,8,10")
    parser.add_argument("--constraints-list", default="3,5,7")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=28)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-26")
    parser.add_argument("--s-order", type=int, default=44)
    parser.add_argument("--u-max", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=130)
    parser.add_argument("--generalized-tol", type=float, default=1e-12)
    parser.add_argument("--close-tol", type=float, default=1e-8)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="volterra_feature_contraction.json")
    args = parser.parse_args()

    bases = parse_ints(args.basis_list)
    constraints = parse_ints(args.constraints_list)
    if len(bases) != len(constraints):
        raise SystemExit("--basis-list and --constraints-list must have the same length")

    rows = [run_case(args, basis, con) for basis, con in zip(bases, constraints)]
    worst = max(rows, key=lambda row: row["topContractionEigenvalue"])
    all_closed = all(row["finiteDominationClosed"] for row in rows)
    data = {
        "theoremName": "Signed Volterra feature contraction",
        "targetInequality": "minus-feature Gram <= plus-feature Gram on Green-minimizer trace image",
        "contractionForm": "lambda_max(P^+ M) <= 1",
        "statuses": {
            "finiteContractionStatus": status(
                "finite signed-feature contraction",
                all_closed,
                (
                    "Every scanned direct-Volterra finite quotient has top "
                    "generalized eigenvalue lambda_max(P^+M) at or below one "
                    "within the displayed tolerance."
                ),
                blocker=None if all_closed else "A finite section has lambda_max(P^+M)>1.",
            ),
            "continuumContractionStatus": status(
                "continuum constrained Volterra contraction",
                False,
                (
                    "The finite spectrum identifies the right inequality, but "
                    "the continuum Hardy/Volterra proof is still needed."
                ),
                blocker=(
                    "Derive the transport/Hardy identity showing the minus "
                    "feature is a contraction of the plus feature on the "
                    "Green-minimizer trace image."
                ),
            ),
        },
        "rows": rows,
        "worstRow": worst,
        "correctedConclusion": (
            "The positive-square problem is exactly the contraction "
            "M <= P for the minus and plus signed Volterra feature Grams. "
            "The finite direct-Volterra quotient scans support the domination "
            "and expose the top trace direction for the analytic Hardy/Green "
            "proof."
        ),
        "nextProofTarget": (
            "Use the top contraction direction to derive a u-transport Hardy "
            "identity for G_minus as a Volterra transform of G_plus, with "
            "boundary terms killed by the Euler-Lagrange Green-minimizer "
            "trace equations."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Signed Volterra feature contraction")
    print(f"  rows: {len(rows)}")
    print(f"  all finite rows closed: {all_closed}")
    for row in rows:
        print(
            f"  basis={row['basis']} constraints={row['constraints']} "
            f"rank={row['traceRank']} top={row['topContractionEigenvalue']:.12e} "
            f"gap={row['topContractionGapToOne']:.12e} "
            f"defect_min={row['defectMin']:.12e}"
        )
    print(f"  worst top contraction: {worst['topContractionEigenvalue']:.12e}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
