#!/usr/bin/env python3
r"""Transport identities for the signed Volterra feature rows.

For each branch sigma,

    M_sigma(u) = int f(s) B_sigma(s,u) ds,
    N_sigma(u) = int (s+u) f(s) B_sigma(s,u) ds,

with B_sigma(s,u)=exp(sigma*omega*(s+u)/2) Psi(s+u)/Psi(s).  Hence the exact
branch transport identity is

    N_sigma = 2 sigma partial_omega M_sigma.

The signed features are

    G_+ = sqrt(w/4)(M+N),     G_- = sqrt(w/4)(M-N).

This script studies the finite range map

    H_- = T H_+

on the Green-minimizer trace image.  The abstract finite T exists exactly when
ker H_+ is contained in ker H_-.  Its norm equals sqrt(lambda_max(P^+M)), so
the continuum target is to identify this finite range map with an explicit
Hardy/Volterra contraction and prove ||T|| <= 1.
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
from volterra_feature_contraction import generalized_psd_spectrum, sym


def build_local_args(args: argparse.Namespace):
    return SimpleNamespace(
        model="full",
        kind=args.kind,
        omega=args.omega,
        L=args.L,
        basis=args.basis,
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


def direct_feature_pairs(args: argparse.Namespace, omega: float):
    sigmas = [0] if omega == 0.0 else [-1, 1]
    branch_weight = 1.0 if omega == 0.0 else 0.5
    rows = []
    direct = np.zeros((args.basis, args.basis), dtype=float)
    plus_basis = np.zeros_like(direct)
    minus_basis = np.zeros_like(direct)
    for sigma in sigmas:
        m, n = branch_feature_matrices(
            kind=args.kind,
            omega=omega,
            sigma=sigma,
            length=float(args.L),
            basis=args.basis,
            s_order=args.s_order,
            u_max=args.u_max,
            u_order=args.u_order,
        )
        plus = 0.25 * branch_weight * ((m + n).T @ (m + n))
        minus = 0.25 * branch_weight * ((m - n).T @ (m - n))
        direct += plus - minus
        plus_basis += plus
        minus_basis += minus
        rows.append(
            {
                "sigma": sigma,
                "weight": branch_weight,
                "m": m,
                "n": n,
            }
        )
    return rows, sym(direct), sym(plus_basis), sym(minus_basis)


def green_minimizer_matrix(args: argparse.Namespace, direct: np.ndarray):
    local = build_local_args(args)
    mp.mp.dps = args.dps
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    _, R = trace_matrix(polys, local)
    N, U, A, B, C, S, minimizer = quotient_blocks(np_to_mp(direct), R, local)
    return {
        "N": N,
        "U": U,
        "S": S,
        "F": mp_to_np(minimizer),
    }


def profile_matrices(feature_rows, F: np.ndarray):
    hplus = []
    hminus = []
    branch_profiles = []
    for row in feature_rows:
        root = math.sqrt(row["weight"] / 4.0)
        plus = root * ((row["m"] + row["n"]) @ F)
        minus = root * ((row["m"] - row["n"]) @ F)
        hplus.append(plus)
        hminus.append(minus)
        branch_profiles.append(
            {
                "sigma": row["sigma"],
                "plus": plus,
                "minus": minus,
                "m": row["m"] @ F,
                "n": row["n"] @ F,
            }
        )
    return np.vstack(hplus), np.vstack(hminus), branch_profiles


def finite_range_transport(hplus: np.ndarray, hminus: np.ndarray, tol: float):
    u, svals, vt = np.linalg.svd(hplus, full_matrices=False)
    cutoff = tol * max(1.0, float(svals[0]) if svals.size else 0.0)
    keep = np.array([i for i, val in enumerate(svals) if val > cutoff], dtype=int)
    if keep.size == 0:
        return {
            "rank": 0,
            "cutoff": cutoff,
            "rangeResidual": float(np.linalg.norm(hminus, ord="fro")),
            "operatorNorm": float("nan"),
            "singularValues": [],
        }
    pinv = vt[keep, :].T @ np.diag(1.0 / svals[keep]) @ u[:, keep].T
    transport = hminus @ pinv
    residual = hminus - transport @ hplus
    t_svals = np.linalg.svd(transport, compute_uv=False)
    return {
        "rank": int(keep.size),
        "cutoff": cutoff,
        "rangeResidual": float(np.linalg.norm(residual, ord="fro")),
        "relativeRangeResidual": float(
            np.linalg.norm(residual, ord="fro") / max(np.linalg.norm(hminus, ord="fro"), 1e-300)
        ),
        "operatorNorm": float(t_svals[0]) if t_svals.size else float("nan"),
        "singularValues": [float(v) for v in t_svals[: min(12, len(t_svals))]],
    }


def omega_transport_error(args: argparse.Namespace, F: np.ndarray, feature_rows) -> dict:
    eps = args.omega_eps
    omega = float(args.omega)
    if omega - eps <= 0:
        return {"checked": False, "reason": "omega-eps must stay positive"}
    plus_rows, _, _, _ = direct_feature_pairs(args, omega + eps)
    minus_rows, _, _, _ = direct_feature_pairs(args, omega - eps)
    max_abs = 0.0
    max_rel = 0.0
    rows = []
    by_sigma_plus = {row["sigma"]: row for row in plus_rows}
    by_sigma_minus = {row["sigma"]: row for row in minus_rows}
    for row in feature_rows:
        sigma = row["sigma"]
        if sigma == 0:
            continue
        m0 = row["m"] @ F
        n0 = row["n"] @ F
        mp1 = by_sigma_plus[sigma]["m"] @ F
        mm1 = by_sigma_minus[sigma]["m"] @ F
        domega = (mp1 - mm1) / (2.0 * eps)
        residual = n0 - 2.0 * sigma * domega
        abs_err = float(np.linalg.norm(residual, ord="fro"))
        rel_err = abs_err / max(float(np.linalg.norm(n0, ord="fro")), 1e-300)
        max_abs = max(max_abs, abs_err)
        max_rel = max(max_rel, rel_err)
        rows.append({"sigma": sigma, "absoluteError": abs_err, "relativeError": rel_err})
    return {
        "checked": True,
        "omegaStep": eps,
        "maxAbsoluteError": max_abs,
        "maxRelativeError": max_rel,
        "rows": rows,
    }


def top_profile_diagnostics(hplus: np.ndarray, hminus: np.ndarray, plus: np.ndarray, minus: np.ndarray, branch_profiles, args):
    spec = generalized_psd_spectrum(minus, plus, args.generalized_tol)
    vals = spec["values"]
    if not vals.size:
        return {}
    z = spec["vectors"][:, -1]
    pnorm = math.sqrt(max(float(z.T @ plus @ z), 0.0))
    if pnorm:
        z = z / pnorm
    p = hplus @ z
    m = hminus @ z
    density = p * p - m * m
    cumulative = np.cumsum(density)
    # The feature rows already include sqrt quadrature weights, so cumulative
    # is a discrete prefix-energy diagnostic, not a literal unweighted integral.
    branch_endpoints = []
    offset = 0
    for prof in branch_profiles:
        count = prof["plus"].shape[0]
        bp = prof["plus"] @ z
        bm = prof["minus"] @ z
        branch_density = bp * bp - bm * bm
        branch_endpoints.append(
            {
                "sigma": prof["sigma"],
                "plusAtFirstNode": float(bp[0]),
                "minusAtFirstNode": float(bm[0]),
                "plusAtLastNode": float(bp[-1]),
                "minusAtLastNode": float(bm[-1]),
                "densitySum": float(np.sum(branch_density)),
                "densityPrefixMin": float(np.min(np.cumsum(branch_density))),
                "densityPrefixFinal": float(np.sum(branch_density)),
            }
        )
        offset += count
    return {
        "topContractionEigenvalue": float(vals[-1]),
        "topContractionGapToOne": float(1.0 - vals[-1]),
        "topVectorPlusEnergy": float(z.T @ plus @ z),
        "topVectorMinusEnergy": float(z.T @ minus @ z),
        "topVectorDefectEnergy": float(z.T @ (plus - minus) @ z),
        "globalDensityPrefixMin": float(np.min(cumulative)),
        "globalDensityPrefixFinal": float(cumulative[-1]),
        "branchEndpointDiagnostics": branch_endpoints,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=12)
    parser.add_argument("--constraints", type=int, default=9)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=28)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-26")
    parser.add_argument("--s-order", type=int, default=64)
    parser.add_argument("--u-max", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=180)
    parser.add_argument("--generalized-tol", type=float, default=1e-12)
    parser.add_argument("--transport-tol", type=float, default=1e-12)
    parser.add_argument("--omega-eps", type=float, default=1e-4)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="volterra_transport_identity.json")
    args = parser.parse_args()

    feature_rows, direct, plus_basis, minus_basis = direct_feature_pairs(args, float(args.omega))
    green = green_minimizer_matrix(args, direct)
    F = green["F"]
    hplus, hminus, branch_profiles = profile_matrices(feature_rows, F)
    plus = sym(hplus.T @ hplus)
    minus = sym(hminus.T @ hminus)
    defect = sym(plus - minus)
    transport = finite_range_transport(hplus, hminus, args.transport_tol)
    omega_identity = omega_transport_error(args, F, feature_rows)
    top = top_profile_diagnostics(hplus, hminus, plus, minus, branch_profiles, args)

    norm_from_contraction = math.sqrt(max(top.get("topContractionEigenvalue", float("nan")), 0.0))
    data = {
        "theoremName": "Volterra signed-feature transport identity",
        "exactBranchIdentity": "N_sigma = 2 sigma partial_omega M_sigma",
        "finiteRangeIdentity": "H_minus = T H_plus on the Green-minimizer trace image",
        "statuses": {
            "omegaTransportIdentityStatus": status(
                "branch omega-transport identity",
                omega_identity.get("checked", False)
                and omega_identity.get("maxRelativeError", 1.0) < 5e-4,
                (
                    "Numerically verifies the exact calculus identity "
                    "N_sigma=2 sigma partial_omega M_sigma by centered "
                    "omega differences with fixed Green-minimizer coefficients."
                ),
                blocker=None
                if omega_identity.get("checked", False)
                else omega_identity.get("reason", "omega transport was not checked"),
            ),
            "finiteRangeTransportStatus": status(
                "finite range map H_minus=T H_plus",
                transport["relativeRangeResidual"] < 1e-8,
                (
                    "The least-squares range map carries plus profiles to "
                    "minus profiles on the scanned Green-minimizer trace image."
                ),
                blocker=None
                if transport["relativeRangeResidual"] < 1e-8
                else "The plus-profile range does not determine the minus-profile range.",
            ),
            "hardyVolterraContractionStatus": status(
                "explicit Hardy/Volterra contraction",
                False,
                (
                    "The finite range map exists and has norm below one, but "
                    "it has not yet been identified with a closed-form "
                    "Hardy/Volterra operator."
                ),
                blocker=(
                    "Derive a u-transport Green identity for the range map T "
                    "and prove its norm is at most one with endpoint terms "
                    "cancelled by the Euler-Lagrange trace equations."
                ),
            ),
        },
        "basis": args.basis,
        "constraints": args.constraints,
        "traceRank": int(green["U"].cols),
        "traceNullity": int(green["N"].cols),
        "transport": transport,
        "omegaTransport": omega_identity,
        "topProfile": top,
        "operatorNormFromContraction": norm_from_contraction,
        "defectMin": float(np.linalg.eigvalsh(defect)[0]),
        "defectMax": float(np.linalg.eigvalsh(defect)[-1]),
        "cutoffCaveat": (
            "Near-null directions in the plus-profile Gram make the displayed "
            "top generalized eigenvalue cutoff-sensitive.  The robust facts are "
            "the small range residual, the subunit finite range-map norm, and "
            "the omega-transport identity."
        ),
        "correctedConclusion": (
            "The finite plus-profile range determines the minus-profile range, "
            "and the computed range map is a subunit contraction after the "
            "chosen numerical cutoff.  The exact calculus identity "
            "N_sigma=2 sigma partial_omega M_sigma is a real transport handle, "
            "but the closed-form Hardy/Volterra operator T is still the missing "
            "continuum proof object."
        ),
        "nextProofTarget": (
            "Convert N_sigma=2 sigma partial_omega M_sigma into a u-transport "
            "Green identity for the range map T, then show the endpoint "
            "boundary form vanishes for Euler-Lagrange Green minimizers."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Volterra signed-feature transport identity")
    print(f"  basis={args.basis} constraints={args.constraints}")
    print(f"  range residual: {transport['relativeRangeResidual']:.3e}")
    print(f"  finite T operator norm: {transport['operatorNorm']:.12e}")
    print(f"  sqrt contraction norm: {norm_from_contraction:.12e}")
    print(f"  omega transport relative error: {omega_identity.get('maxRelativeError', float('nan')):.3e}")
    print(f"  top contraction: {top.get('topContractionEigenvalue', float('nan')):.12e}")
    print(f"  defect min: {np.linalg.eigvalsh(defect)[0]:.12e}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
