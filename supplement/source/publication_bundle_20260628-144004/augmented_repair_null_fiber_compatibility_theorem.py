#!/usr/bin/env python3
r"""Null-fiber compatibility for augmented repair forms."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-consequence-json",
        default="xi_augmented_finite_schur_repair_consequence_theorem.json",
    )
    parser.add_argument(
        "--trace-quotient-json",
        default="augmented_trace_quotient_compatibility_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_repair_null_fiber_compatibility_theorem.json",
    )
    args = parser.parse_args()

    finite = load(args.finite_consequence_json)
    quotient = load(args.trace_quotient_json)

    finite_mu_ok = bool(finite.get("finiteMuKernelAnnihilationClosed"))
    finite_range_ok = bool(finite.get("finiteSchurRangeConditionClosed"))
    quotient_ok = bool(quotient.get("traceQuotientCompatibilityClosed"))
    compatibility_ok = finite_mu_ok and finite_range_ok and quotient_ok

    data = {
        "theoremName": "augmented repair null-fiber compatibility theorem",
        "proofClass": "symbolic identity with interval/ball certificate inputs",
        "finiteConsequenceJson": args.finite_consequence_json,
        "traceQuotientJson": args.trace_quotient_json,
        "statement": (
            "The augmented repair forms are compatible with the finite null "
            "fibers and hence with the transported trace quotient relation."
        ),
        "statuses": {
            "finiteMuNullFiberInputStatus": status(
                "finite Mu null-fiber input",
                finite_mu_ok,
                "The finite Schur consequence proves Mu annihilates ker R_aug,N.",
            ),
            "finiteSchurRangeInputStatus": status(
                "finite Schur range input",
                finite_range_ok,
                "The finite Schur range condition makes the repair depend only on augmented trace coordinates.",
            ),
            "traceQuotientCompatibilityInputStatus": status(
                "trace quotient compatibility input",
                quotient_ok,
                "The quotient compatibility theorem identifies finite trace equivalence with X_aug.",
            ),
            "repairNullFiberCompatibilityStatus": status(
                "repair null-fiber compatibility",
                compatibility_ok,
                (
                    "On the finite cores the repair vanishes on trace-null "
                    "directions; quotient compatibility transports the same "
                    "representative independence to X_aug."
                ),
            ),
        },
        "repairNullFiberCompatibilityClosed": compatibility_ok,
        "finiteRepairNullFiberCompatibilityClosed": finite_mu_ok and finite_range_ok,
        "transportedRepairRepresentativeIndependenceClosed": compatibility_ok,
        "proof": [
            "Import finite Mu annihilation and Schur range compatibility.",
            "Those finite identities say the repair value is unchanged by adding a vector in ker R_aug,N.",
            "Use trace quotient compatibility to pass representative independence to the transported trace quotient.",
        ],
        "remainingAnalyticGap": None
        if compatibility_ok
        else "Finite null-fiber compatibility or trace quotient compatibility is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented repair null-fiber compatibility theorem")
    print(f"  finite Mu null-fiber: {finite_mu_ok}")
    print(f"  finite Schur range: {finite_range_ok}")
    print(f"  quotient compatibility: {quotient_ok}")
    print(f"  repair compatibility: {compatibility_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
