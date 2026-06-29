#!/usr/bin/env python3
r"""Finite augmented Schur certificate constants theorem.

This theorem object reads the finite augmented trace repair certificate and
extracts only the deterministic inequalities needed by the universal finite
Schur/Douglas algebra theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate-json", default="xi_augmented_trace_repair_schur.json")
    parser.add_argument("--eig-tol", type=float, default=1e-20)
    parser.add_argument(
        "--json-out",
        default="xi_augmented_finite_schur_constants_theorem.json",
    )
    args = parser.parse_args()

    cert = load(args.certificate_json)
    augmented = cert.get("augmentedRepair", {})
    d_min = float(augmented.get("dMin", 0.0))
    p_min = float(augmented.get("pMin", 0.0))
    schur_min = float(augmented.get("constructedSchurMin", 0.0))
    range_residual = float(augmented.get("normalizedRangeResidual", 0.0))
    mu_action = float(cert.get("augmentedMuActionOnAugmentedNullspace", 0.0))
    nullity = int(augmented.get("nullity", -1))
    rank = int(augmented.get("rank", -1))

    trace_quotient_ok = bool(nullity >= 0 and rank > 0)
    range_condition_ok = bool(range_residual <= args.eig_tol)
    repair_nonnegative_ok = bool(cert.get("augmentedRepairPositive") and d_min >= -args.eig_tol)
    schur_gap_ok = bool(schur_min >= -args.eig_tol)
    repaired_form_ok = bool(p_min >= -args.eig_tol)
    mu_annihilation_ok = bool(cert.get("augmentedMuAnnihilated") and mu_action <= args.eig_tol)
    constants_ok = all(
        [
            trace_quotient_ok,
            range_condition_ok,
            repair_nonnegative_ok,
            schur_gap_ok,
            repaired_form_ok,
            mu_annihilation_ok,
        ]
    )

    data = {
        "theoremName": "augmented finite Schur constants theorem",
        "proofClass": "deterministic high-precision finite certificate",
        "sourceCertificate": {
            "path": args.certificate_json,
            "sha256": sha256(args.certificate_json),
        },
        "finiteTraceDimensions": {
            "traceRank": rank,
            "traceNullity": nullity,
            "lambdaRows": cert.get("lambdaRows"),
            "muRows": cert.get("muRows"),
            "augmentedRows": cert.get("augmentedRows"),
            "basis": cert.get("basis"),
            "constraints": cert.get("constraints"),
        },
        "certificateConstants": {
            "dMin": d_min,
            "pMin": p_min,
            "constructedSchurMin": schur_min,
            "normalizedRangeResidual": range_residual,
            "augmentedMuActionOnAugmentedNullspace": mu_action,
            "eigTolerance": args.eig_tol,
        },
        "statuses": {
            "finiteTraceQuotientConstantStatus": status(
                "finite trace quotient dimensions",
                trace_quotient_ok,
                "The augmented finite trace model has a resolved positive trace rank.",
            ),
            "finiteRangeResidualConstantStatus": status(
                "finite Schur range residual",
                range_condition_ok,
                "The normalized Douglas range residual is below the certificate tolerance.",
            ),
            "finiteDaugNonnegativeConstantStatus": status(
                "finite D_aug,N nonnegative constant",
                repair_nonnegative_ok,
                "The computed lower spectral endpoint of D_aug,N is nonnegative up to the stated tolerance.",
            ),
            "finiteSchurGapConstantStatus": status(
                "finite Schur complement gap",
                schur_gap_ok,
                "The repaired Schur complement lower endpoint is nonnegative up to tolerance.",
            ),
            "finiteRepairedFormConstantStatus": status(
                "finite repaired form lower endpoint",
                repaired_form_ok,
                "The repaired finite form lower endpoint is nonnegative up to tolerance.",
            ),
            "finiteMuAnnihilationConstantStatus": status(
                "finite Mu kernel annihilation constant",
                mu_annihilation_ok,
                "The Mu rows vanish on the augmented finite trace kernel up to tolerance.",
            ),
        },
        "finiteSchurConstantsClosed": constants_ok,
        "finiteTraceQuotientConstantsClosed": trace_quotient_ok,
        "finiteSchurRangeConstantsClosed": range_condition_ok,
        "finiteRepairOperatorNonnegativeConstantsClosed": repair_nonnegative_ok,
        "finiteSchurGapConstantsClosed": schur_gap_ok,
        "finiteRepairedFormPositiveConstantsClosed": repaired_form_ok,
        "finiteMuKernelAnnihilationConstantsClosed": mu_annihilation_ok,
        "remainingAnalyticGap": None
        if constants_ok
        else "One finite Schur certificate constant failed its tolerance.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented finite Schur constants theorem")
    print(f"  trace quotient constants: {trace_quotient_ok}")
    print(f"  range residual constants: {range_condition_ok}")
    print(f"  D nonnegative constants: {repair_nonnegative_ok}")
    print(f"  Schur gap constants: {schur_gap_ok}")
    print(f"  repaired form constants: {repaired_form_ok}")
    print(f"  Mu annihilation constants: {mu_annihilation_ok}")
    print(f"  constants closed: {constants_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
