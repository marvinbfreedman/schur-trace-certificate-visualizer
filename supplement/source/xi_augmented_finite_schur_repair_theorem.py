#!/usr/bin/env python3
r"""Legacy wrapper for the finite augmented Schur repair theorem.

The proof has been split into the universal finite Schur/Douglas algebra
theorem and the finite augmented constants theorem.  The narrow consequence
file used upstream is
``xi_augmented_finite_schur_repair_consequence_theorem.json``.
"""

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
        "--consequence-json",
        default="xi_augmented_finite_schur_repair_consequence_theorem.json",
    )
    parser.add_argument("--certificate-json", default="")
    parser.add_argument("--eig-tol", type=float, default=1e-20)
    parser.add_argument("--json-out", default="xi_augmented_finite_schur_repair_theorem.json")
    args = parser.parse_args()

    consequence = load(args.consequence_json)
    repair_nonnegative = bool(consequence.get("finiteRepairOperatorNonnegativeClosed"))
    repaired_form_positive = bool(consequence.get("finiteRepairedFormPositiveClosed"))
    range_condition = bool(consequence.get("finiteSchurRangeConditionClosed"))
    mu_annihilation = bool(consequence.get("finiteMuKernelAnnihilationClosed"))
    trace_quotient_ok = bool(
        (consequence.get("finiteTraceDimensions") or {}).get("traceRank", 0) > 0
    )
    finite_closed = bool(consequence.get("finiteAugmentedSchurRepairClosed"))

    data = {
        "theoremName": "finite augmented Schur repair theorem legacy wrapper",
        "proofClass": "symbolic finite consequence wrapper",
        "finiteSchurRepairConsequenceJson": args.consequence_json,
        "legacyInputs": {
            "certificateJson": args.certificate_json or None,
            "eigTolerance": args.eig_tol,
        },
        "finiteTraceDimensions": consequence.get("finiteTraceDimensions"),
        "certificateConstants": consequence.get("certificateConstants"),
        "statuses": {
            "finiteTraceQuotientStatus": status(
                "finite augmented trace quotient",
                trace_quotient_ok,
                (
                    "The augmented trace map has a resolved trace image and a "
                    "well-defined finite quotient/kernel splitting."
                ),
            ),
            "rangeConditionStatus": status(
                "finite Schur range condition",
                range_condition,
                (
                    "The normalized range residual for the off-diagonal Schur "
                    "block is below tolerance."
                ),
            ),
            "repairOperatorNonnegativeStatus": status(
                "finite repair operator nonnegative",
                repair_nonnegative,
                "The constructed D_aug,N has no resolved negative spectrum.",
            ),
            "repairedFormPositiveStatus": status(
                "finite repaired form positive",
                repaired_form_positive,
                (
                    "The constructed Schur complement and the repaired finite "
                    "form K+R_aug^*D_aug,N R_aug are nonnegative."
                ),
            ),
            "muKernelAnnihilationStatus": status(
                "Mu annihilates finite augmented trace kernel",
                mu_annihilation,
                (
                    "The primitive Mu rows vanish on ker R_aug in the finite "
                    "quotient model."
                ),
            ),
        },
        "finiteAugmentedSchurRepairClosed": finite_closed,
        "finiteRepairOperatorNonnegativeClosed": repair_nonnegative,
        "finiteRepairedFormPositiveClosed": repaired_form_positive,
        "finiteMuKernelAnnihilationClosed": mu_annihilation,
        "finiteSchurRangeConditionClosed": range_condition,
        "formalProof": [
            "Import the narrow finite augmented Schur repair consequence theorem.",
            "Re-export legacy compatibility flags expected by older theorem scripts.",
        ],
        "remainingAnalyticGap": None if finite_closed else "One finite Schur repair condition is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Finite augmented Schur repair theorem")
    print(f"  trace quotient: {trace_quotient_ok}")
    print(f"  range condition: {range_condition}")
    print(f"  D nonnegative: {repair_nonnegative}")
    print(f"  repaired form positive: {repaired_form_positive}")
    print(f"  Mu kernel annihilation: {mu_annihilation}")
    print(f"  theorem closed: {finite_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
