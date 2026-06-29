#!/usr/bin/env python3
r"""Hardy multiplier reduction for the signed Volterra feature contraction.

For r=s+u >= 0, the signed features are generated from the same lifted
integrand:

    G_+(u) = int sqrt(w/4) (1+r) f(s) B_sigma(s,u) ds,
    G_-(u) = int sqrt(w/4) (1-r) f(s) B_sigma(s,u) ds.

Thus before the Volterra s-integration,

    lifted_minus = kappa(r) lifted_plus,
    kappa(r) = (1-r)/(1+r),        |kappa(r)| <= 1.

The actual finite range map is therefore a compression

    T = C K E,

where C integrates over s, K is multiplication by kappa, and E is the
Green-minimizer right inverse/lift from the observed plus profile to the lifted
plus integrand.  Pointwise |K|<=1 is exact; the remaining continuum theorem is
that the Euler-Lagrange Green lift E is the correct minimal lift and kills the
endpoint boundary form, making the compression contractive.
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
from reduced_exact_finite import endpoint_ratios, pieces
from trace_volterra_green_feature_map import (
    mp_to_np,
    np_to_mp,
    quadrature,
    quotient_blocks,
    shifted_legendre_values,
    status,
)
from volterra_transport_identity import finite_range_transport


def sym(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat + mat.T)


def a_value(s: float, u: float, pcs) -> float:
    es = math.exp(s)
    eu = math.exp(u)
    total = 0.0
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * math.exp(beta * u - c * es * (eu - 1.0))
    return total


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


def lifted_feature_matrices(args: argparse.Namespace):
    pcs = pieces(args.kind)
    sigmas = [0] if float(args.omega) == 0.0 else [-1, 1]
    branch_weight = 1.0 if float(args.omega) == 0.0 else 0.5
    s_pts, s_wts = quadrature(float(args.L), args.s_order)
    u_pts, u_wts = quadrature(args.u_max, args.u_order)
    phi = shifted_legendre_values(args.basis, float(args.L), s_pts)
    rows_plus = []
    rows_minus = []
    rows_kplus = []
    kappas = []
    endpoint_rows = []
    for sigma in sigmas:
        root_branch = math.sqrt(branch_weight / 4.0)
        first_plus = np.zeros(args.basis, dtype=float)
        first_minus = np.zeros(args.basis, dtype=float)
        last_plus = np.zeros(args.basis, dtype=float)
        last_minus = np.zeros(args.basis, dtype=float)
        for iu, u in enumerate(u_pts):
            plus_integrated = np.zeros(args.basis, dtype=float)
            minus_integrated = np.zeros(args.basis, dtype=float)
            kplus_integrated = np.zeros(args.basis, dtype=float)
            for js, s in enumerate(s_pts):
                r = float(s + u)
                branch = math.exp(0.5 * sigma * float(args.omega) * r)
                base = (
                    root_branch
                    * math.sqrt(float(u_wts[iu]))
                    * float(s_wts[js])
                    * a_value(float(s), float(u), pcs)
                    * branch
                )
                kappa = (1.0 - r) / (1.0 + r)
                plus = base * (1.0 + r) * phi[js, :]
                minus = base * (1.0 - r) * phi[js, :]
                plus_integrated += plus
                minus_integrated += minus
                kplus_integrated += kappa * plus
                kappas.append(kappa)
            rows_plus.append(plus_integrated)
            rows_minus.append(minus_integrated)
            rows_kplus.append(kplus_integrated)
            if iu == 0:
                first_plus = plus_integrated
                first_minus = minus_integrated
            if iu == len(u_pts) - 1:
                last_plus = plus_integrated
                last_minus = minus_integrated
        endpoint_rows.append(
            {
                "sigma": sigma,
                "firstU": float(u_pts[0]),
                "lastU": float(u_pts[-1]),
                "_firstPlus": first_plus,
                "_firstMinus": first_minus,
                "_lastPlus": last_plus,
                "_lastMinus": last_minus,
            }
        )
    return (
        np.vstack(rows_plus),
        np.vstack(rows_minus),
        np.vstack(rows_kplus),
        np.array(kappas, dtype=float),
        endpoint_rows,
    )


def green_minimizer(args: argparse.Namespace, plus_basis: np.ndarray, minus_basis: np.ndarray):
    direct = sym(plus_basis - minus_basis)
    local = build_local_args(args)
    mp.mp.dps = args.dps
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    _, R = trace_matrix(polys, local)
    N, U, A, B, C, S, minimizer = quotient_blocks(np_to_mp(direct), R, local)
    return mp_to_np(minimizer), int(U.cols), int(N.cols)


def endpoint_pairings(endpoint_rows, F: np.ndarray):
    out = []
    for row in endpoint_rows:
        fp = row["_firstPlus"] @ F
        fm = row["_firstMinus"] @ F
        lp = row["_lastPlus"] @ F
        lm = row["_lastMinus"] @ F
        out.append(
            {
                "sigma": row["sigma"],
                "firstU": row["firstU"],
                "lastU": row["lastU"],
                "firstPlusNorm": float(np.linalg.norm(fp)),
                "firstMinusNorm": float(np.linalg.norm(fm)),
                "lastPlusNorm": float(np.linalg.norm(lp)),
                "lastMinusNorm": float(np.linalg.norm(lm)),
                "firstPairingTrace": float(np.dot(fp, fm)),
                "lastPairingTrace": float(np.dot(lp, lm)),
            }
        )
    return out


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
    parser.add_argument("--transport-tol", type=float, default=1e-12)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="volterra_hardy_transport_derivation.json")
    args = parser.parse_args()

    plus_basis_rows, minus_basis_rows, kplus_rows, kappas, endpoint_rows = lifted_feature_matrices(args)
    plus_basis = sym(plus_basis_rows.T @ plus_basis_rows)
    minus_basis = sym(minus_basis_rows.T @ minus_basis_rows)
    F, trace_rank, trace_nullity = green_minimizer(args, plus_basis, minus_basis)
    hplus = plus_basis_rows @ F
    hminus = minus_basis_rows @ F
    hkplus = kplus_rows @ F
    transport = finite_range_transport(hplus, hminus, args.transport_tol)
    k_identity_error = np.linalg.norm(hminus - hkplus, ord="fro")
    k_identity_relative = k_identity_error / max(np.linalg.norm(hminus, ord="fro"), 1e-300)
    endpoint = endpoint_pairings(endpoint_rows, F)
    kappa_max = float(np.max(np.abs(kappas)))
    kappa_min = float(np.min(kappas))
    kappa_max_raw = float(np.max(kappas))

    data = {
        "theoremName": "Hardy multiplier reduction for signed Volterra transport",
        "explicitMultiplier": "kappa(s,u)=(1-s-u)/(1+s+u)",
        "compressedOperator": "T = C K E on the Green-minimizer plus-profile range",
        "statuses": {
            "pointwiseHardyMultiplierStatus": status(
                "pointwise Hardy multiplier",
                bool(kappa_max <= 1.0 + 1e-15),
                "For s,u>=0, |(1-s-u)/(1+s+u)|<=1 exactly.",
            ),
            "compressedIdentityStatus": status(
                "compressed multiplier identity",
                bool(k_identity_relative < 1e-10),
                "The integrated minus profiles equal C K applied to the lifted plus profiles.",
                blocker=None
                if k_identity_relative < 1e-10
                else "The lifted multiplier identity failed numerically.",
            ),
            "greenMinimalLiftStatus": status(
                "Green-minimizer minimal lift",
                False,
                (
                    "This is the remaining analytic theorem: the Euler-Lagrange "
                    "Green lift must be shown to be the minimal contraction lift "
                    "for the compressed Hardy multiplier."
                ),
                blocker=(
                    "Use the trace-fiber Euler-Lagrange equation to show the "
                    "boundary concomitant from integrating C K E by parts "
                    "vanishes on the completed trace image."
                ),
            ),
        },
        "basis": args.basis,
        "constraints": args.constraints,
        "traceRank": trace_rank,
        "traceNullity": trace_nullity,
        "kappaMin": kappa_min,
        "kappaMax": kappa_max_raw,
        "kappaAbsMax": kappa_max,
        "compressedIdentityFrobenius": float(k_identity_error),
        "compressedIdentityRelative": float(k_identity_relative),
        "transport": transport,
        "endpointDiagnostics": endpoint,
        "correctedConclusion": (
            "The explicit pointwise Hardy multiplier is identified and exactly "
            "reproduces the integrated minus features from the lifted plus "
            "features.  The remaining step is not algebraic identification of "
            "T; it is proving that the Green-minimizer right inverse E makes "
            "the compression C K E contractive, with boundary terms killed by "
            "the Euler-Lagrange trace equations."
        ),
        "nextProofTarget": (
            "Write the integration-by-parts formula for C K E in the lifted "
            "(s,u) variables and express its boundary concomitant as the "
            "Euler-Lagrange trace-fiber residual Q(f_x,h), which is zero for "
            "all h in ker R."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Hardy multiplier reduction")
    print(f"  kappa abs max: {kappa_max:.12e}")
    print(f"  compressed identity rel error: {k_identity_relative:.3e}")
    print(f"  finite T operator norm: {transport['operatorNorm']:.12e}")
    print(f"  range residual: {transport['relativeRangeResidual']:.3e}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
