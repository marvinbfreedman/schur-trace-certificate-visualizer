#!/usr/bin/env python3
r"""Scan the quotient repair-vanishing condition.

After the primitive boundary audit, the original lift can no longer rely on a
positive endpoint repair: the primitive boundary operator is D_bdy=0.  The
remaining repair route is therefore

    D_q = 0 on X_R,

equivalently

    Gamma^* Gamma - C <= 0

in the quotient Schur coordinates.  This script computes the finite Galerkin
analogue of that operator inequality by reusing quotient_factorization_mp.

The output is finite evidence only.  A negative finite top eigenvalue of
H=Gamma^*Gamma-C suggests the repair-free theorem is plausible; a positive
stable top eigenvalue would disprove this route.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp

from quotient_factorization_mp import gram_matrix, trace_matrix, quotient_certificate


def as_float(x) -> float:
    return float(x)


def fmt(x, digits: int = 12) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def parse_words(text: str) -> list[str]:
    return [piece for piece in text.replace(",", " ").split() if piece]


def build_args(base: argparse.Namespace, *, model: str, basis: int, constraints: int):
    return SimpleNamespace(
        model=model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=max(base.min_quad, base.quad_factor * basis),
        laguerre=base.laguerre,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
        zero_tol=base.zero_tol,
        dps=base.dps,
    )


def scan_one(local) -> dict:
    omega = mp.mpf(local.omega)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, omega, length)
    centers, R = trace_matrix(polys, local)
    cert = quotient_certificate(K, R, local)
    h_max = cert["h_max"]
    h_min = cert["h_min"]
    dq_positive = max(mp.mpf("0"), h_max)
    return {
        "model": local.model,
        "kind": local.kind,
        "omega": as_float(omega),
        "L": as_float(length),
        "basis": local.basis,
        "quad": local.quad,
        "laguerre": local.laguerre,
        "constraints": local.constraints,
        "traceRank": int(cert["rank"]),
        "traceNullity": int(cert["nullity"]),
        "kMin": as_float(cert["kmin"]),
        "aMin": as_float(cert["amin"]),
        "bNorm": as_float(cert["b_norm"]),
        "rangeResidual": as_float(cert["range_resid"]),
        "normalizedRangeResidual": as_float(cert["normalized_range_resid"]),
        "gamma2Min": as_float(cert["gamma2_min"]),
        "gamma2Max": as_float(cert["gamma2_max"]),
        "schurDefectMin": as_float(h_min),
        "schurDefectMax": as_float(h_max),
        "dqPositivePartNorm": as_float(dq_positive),
        "finiteDqZero": h_max <= mp.mpf(local.zero_tol),
        "lowSchurDefect": [as_float(v) for v in cert["low_H"]],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", default="full")
    parser.add_argument("--kind", default="tilde3", choices=["raw1", "raw2", "raw3", "tilde3"])
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis-values", default="6 8")
    parser.add_argument("--constraint-values", default="3 5")
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
    parser.add_argument("--margin", default="0")
    parser.add_argument("--zero-tol", default="1e-20")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="dq_vanishing_schur_defect_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    models = parse_words(args.models)
    basis_values = parse_ints(args.basis_values)
    constraint_values = parse_ints(args.constraint_values)

    rows = []
    print("D_q vanishing Schur-defect scan", flush=True)
    print(
        "  model kind basis constraints rank nullity "
        "h_min h_max dq_pos finite_zero",
        flush=True,
    )
    for model in models:
        for basis in basis_values:
            for constraints in constraint_values:
                if constraints >= basis:
                    continue
                local = build_args(args, model=model, basis=basis, constraints=constraints)
                row = scan_one(local)
                rows.append(row)
                print(
                    f"  {model:>6} {args.kind:>6} {basis:5d} {constraints:11d} "
                    f"{row['traceRank']:4d} {row['traceNullity']:7d} "
                    f"{fmt(row['schurDefectMin'], 8):>12} "
                    f"{fmt(row['schurDefectMax'], 8):>12} "
                    f"{fmt(row['dqPositivePartNorm'], 8):>12} "
                    f"{row['finiteDqZero']}",
                    flush=True,
                )

    worst = max(rows, key=lambda row: row["schurDefectMax"]) if rows else None
    all_zero = bool(rows) and all(row["finiteDqZero"] for row in rows)
    data = {
        "theoremName": "finite D_q vanishing Schur-defect scan",
        "meaning": (
            "Finite Galerkin test of D_q=0, equivalent to "
            "Gamma^*Gamma-C<=0 in quotient Schur coordinates."
        ),
        "parameters": {
            "models": models,
            "kind": args.kind,
            "omega": float(mp.mpf(args.omega)),
            "L": float(mp.mpf(args.L)),
            "basisValues": basis_values,
            "constraintValues": constraint_values,
            "zeroTol": float(mp.mpf(args.zero_tol)),
            "dps": args.dps,
        },
        "rows": rows,
        "worstRow": worst,
        "allFiniteDqZero": all_zero,
        "finiteEvidenceStatus": {
            "label": "finite Schur defect nonpositive",
            "closed": all_zero,
            "status": "closed" if all_zero else "open",
            "reason": (
                "Every scanned finite quotient has lambda_max(Gamma^*Gamma-C) "
                "below the zero tolerance."
                if all_zero
                else "At least one scanned finite quotient has positive Schur defect."
            ),
        },
        "continuumTheoremClosed": False,
        "nextProofTarget": (
            "Promote the finite inequality Gamma^*Gamma-C<=0 to the completed "
            "trace range X_R, or identify a stable positive Schur-defect mode "
            "and abandon the repair-free route."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    if worst:
        print(
            f"  worst h_max={fmt(worst['schurDefectMax'], 12)} "
            f"at model={worst['model']} basis={worst['basis']} "
            f"constraints={worst['constraints']}",
            flush=True,
        )
    print(f"  all finite D_q zero: {all_zero}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
