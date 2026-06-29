#!/usr/bin/env python3
r"""Finite Green-lift boundary/minimality theorem.

The explicit Hardy multiplier reduction gives

    D_trace = P - M,

where P and M are the plus/minus Volterra feature Grams.  For a fixed trace
vector x, the Green lift f_x is the minimizer of Q(f)=<f,(P-M)f> in the fiber
Rf=x.  Therefore every admissible variation h in ker R satisfies

    Q(f_x,h)=0.

Written in feature variables this is the boundary/concomitant cancellation

    <G_+(f_x), G_+(h)> - <G_-(f_x), G_-(h)> = 0.

This is the finite Galerkin form of the continuum target: the boundary
concomitant produced by the compressed Hardy multiplier C K E is exactly the
Euler-Lagrange trace-fiber residual.  The remaining non-finite work is to pass
this identity through the continuum closure and the lifted integration by
parts.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp
import numpy as np

from quotient_factorization_mp import shifted_legendre_polys, trace_matrix
from trace_volterra_green_feature_map import mp_to_np, np_to_mp, quotient_blocks, status
from volterra_hardy_transport_derivation import lifted_feature_matrices, sym


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


def deterministic_vectors(cols: int, count: int) -> list[np.ndarray]:
    out = []
    for seed in range(count):
        v = np.zeros((cols, 1), dtype=float)
        for i in range(cols):
            value = ((seed + 3) * (i + 1)) / 17.0
            if (seed + i) % 2:
                value *= -1.0
            v[i, 0] = value
        out.append(v)
    return out


def scalar(value: np.ndarray) -> float:
    return float(np.asarray(value).reshape(()))


def run_case(args: argparse.Namespace) -> dict:
    mp.mp.dps = args.dps
    plus_rows, minus_rows, _kplus_rows, _kappas, _endpoint_rows = lifted_feature_matrices(args)
    plus = sym(plus_rows.T @ plus_rows)
    minus = sym(minus_rows.T @ minus_rows)
    qmat = sym(plus - minus)

    local = build_local_args(args)
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    _, R = trace_matrix(polys, local)
    N, U, A, B, C, S, minimizer = quotient_blocks(np_to_mp(qmat), R, local)

    Nn = mp_to_np(N)
    Un = mp_to_np(U)
    Fn = mp_to_np(minimizer)
    Sn = mp_to_np(S)
    Rn = mp_to_np(R)

    # Boundary/concomitant matrix: all trace Green lifts against all ker R
    # variations.  This equals F^T(P-M)N.
    plus_cross = Fn.T @ plus @ Nn
    minus_cross = Fn.T @ minus @ Nn
    concomitant = plus_cross - minus_cross
    q_cross = Fn.T @ qmat @ Nn
    el_block = mp_to_np(A) @ (-mp_to_np(B))  # not used for residual; kept structural.

    # Since minimizer=U-N A^+B, the direct residual is F^T Q N.
    concomitant_norm = float(np.linalg.norm(concomitant, ord="fro"))
    q_cross_norm = float(np.linalg.norm(q_cross, ord="fro"))
    plus_cross_norm = float(np.linalg.norm(plus_cross, ord="fro"))
    minus_cross_norm = float(np.linalg.norm(minus_cross, ord="fro"))
    relative = concomitant_norm / max(plus_cross_norm + minus_cross_norm, 1e-300)

    trace_error = np.linalg.norm(Rn @ Fn - Rn @ Un, ord="fro")
    rows = []
    min_gap = float("inf")
    max_gap_error = 0.0
    for idx, (z, y) in enumerate(
        zip(deterministic_vectors(Fn.shape[1], args.test_vectors), deterministic_vectors(Nn.shape[1], args.test_vectors))
    ):
        f = Fn @ z
        h = Nn @ y
        plus_f = plus_rows @ f
        minus_f = minus_rows @ f
        plus_h = plus_rows @ h
        minus_h = minus_rows @ h
        boundary = scalar(plus_f.T @ plus_h - minus_f.T @ minus_h)
        qfh = scalar(f.T @ qmat @ h)
        # Minimality: Q(f+h)-Q(f)=Q(h), because Q(f,h)=0.
        qf = scalar(f.T @ qmat @ f)
        qh = scalar(h.T @ qmat @ h)
        qfh_total = scalar((f + h).T @ qmat @ (f + h))
        gap_error = abs((qfh_total - qf) - qh)
        min_gap = min(min_gap, qfh_total - qf)
        max_gap_error = max(max_gap_error, gap_error)
        rows.append(
            {
                "seed": idx,
                "boundaryConcomitant": boundary,
                "qfh": qfh,
                "minimalityGap": qfh_total - qf,
                "kerEnergy": qh,
                "gapMinusKerEnergy": (qfh_total - qf) - qh,
            }
        )

    svals = np.linalg.eigvalsh(sym(Sn))
    qn_vals = np.linalg.eigvalsh(sym(Nn.T @ qmat @ Nn)) if Nn.size else np.array([0.0])

    return {
        "basis": args.basis,
        "constraints": args.constraints,
        "traceRank": int(U.cols),
        "traceNullity": int(N.cols),
        "schurMin": float(svals[0]),
        "schurMax": float(svals[-1]),
        "kerEnergyMin": float(qn_vals[0]),
        "kerEnergyMax": float(qn_vals[-1]),
        "traceError": float(trace_error),
        "plusCrossNorm": plus_cross_norm,
        "minusCrossNorm": minus_cross_norm,
        "boundaryConcomitantFrobenius": concomitant_norm,
        "qCrossFrobenius": q_cross_norm,
        "relativeBoundaryConcomitant": float(relative),
        "maxMinimalityGapError": float(max_gap_error),
        "minPerturbationGap": float(min_gap),
        "rows": rows,
        "finiteBoundaryTheoremClosed": bool(
            relative < args.close_tol
            and trace_error < args.close_tol
            and max_gap_error < args.energy_tol
            and qn_vals[0] > -args.energy_tol
        ),
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
    parser.add_argument("--test-vectors", type=int, default=5)
    parser.add_argument("--close-tol", type=float, default=1e-8)
    parser.add_argument("--energy-tol", type=float, default=1e-8)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="green_lift_boundary_theorem.json")
    args = parser.parse_args()

    check = run_case(args)
    data = {
        "theoremName": "Green-lift boundary/minimality theorem",
        "finiteIdentity": "<G_+(f_x),G_+(h)>-<G_-(f_x),G_-(h)>=Q(f_x,h)=0 for h in ker R",
        "statuses": {
            "finiteBoundaryConcomitantStatus": status(
                "finite boundary concomitant equals Euler-Lagrange residual",
                check["finiteBoundaryTheoremClosed"],
                (
                    "The plus/minus feature cross-form against ker R equals "
                    "the trace-fiber Euler-Lagrange residual and vanishes for "
                    "the Green minimizer."
                ),
                blocker=None
                if check["finiteBoundaryTheoremClosed"]
                else "Finite boundary residual or perturbation minimality did not close.",
            ),
            "continuumBoundaryTheoremStatus": status(
                "continuum Green-lift boundary theorem",
                False,
                (
                    "The finite Galerkin identity is closed.  The remaining "
                    "work is passing the lifted integration-by-parts "
                    "concomitant to the completed Volterra trace domain."
                ),
                blocker=(
                    "Prove density/closure for the lifted C K E "
                    "integration-by-parts formula and identify its endpoint "
                    "row with Q(f_x,h) for all h in ker R."
                ),
            ),
        },
        "check": check,
        "correctedConclusion": (
            "In finite direct-Volterra Galerkin form, the boundary concomitant "
            "from the plus/minus Hardy split is exactly the Euler-Lagrange "
            "orthogonality Q(f_x,h)=0.  This proves the desired cancellation "
            "at the finite quotient level; the only remaining step is the "
            "continuum closure of the lifted integration-by-parts identity."
        ),
        "nextProofTarget": (
            "Prove the continuum closure theorem: smooth compactly supported "
            "lifts are dense in the completed trace-fiber domain, the C K E "
            "boundary concomitant is continuous in the Volterra form norm, "
            "and its finite identity passes to Q(f_x,h)=0 on ker R."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Green-lift boundary/minimality theorem")
    print(f"  finite boundary theorem closed: {check['finiteBoundaryTheoremClosed']}")
    print(f"  boundary concomitant rel: {check['relativeBoundaryConcomitant']:.3e}")
    print(f"  boundary concomitant frob: {check['boundaryConcomitantFrobenius']:.3e}")
    print(f"  trace error: {check['traceError']:.3e}")
    print(f"  max minimality gap error: {check['maxMinimalityGapError']:.3e}")
    print(f"  ker energy min: {check['kerEnergyMin']:.3e}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
