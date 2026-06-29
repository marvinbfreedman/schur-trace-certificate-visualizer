#!/usr/bin/env python3
r"""Fixed-core finite Schur repaired positivity consequence.

This wrapper exposes only the finite-core input needed by
``augmented_pointwise_repaired_form_nonnegative_certificate.json``:

    q_N^Schur >= 0

for each fixed augmented trace-core index N.

The broader finite Schur repair theorem also exports range, D_aug, and Mu
facts for other layers; those are intentionally not imported here.
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
        "--finite-consequence-json",
        default="xi_augmented_finite_schur_repair_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_fixed_core_schur_positive_consequence_theorem.json",
    )
    args = parser.parse_args()

    finite = load(args.finite_consequence_json)
    repaired_form_ok = bool(finite.get("finiteRepairedFormPositiveClosed"))
    repair_operator_ok = bool(finite.get("finiteRepairOperatorNonnegativeClosed"))
    schur_repair_ok = bool(finite.get("finiteAugmentedSchurRepairClosed"))
    ok = repaired_form_ok and repair_operator_ok and schur_repair_ok

    data = {
        "theoremName": "augmented fixed-core Schur positivity consequence",
        "proofClass": "symbolic identity",
        "finiteConsequenceJson": args.finite_consequence_json,
        "statement": (
            "For each fixed augmented trace-core index N, the finite "
            "Schur-repaired core quadratic form q_N^Schur is nonnegative."
        ),
        "statuses": {
            "finiteSchurRepairInputStatus": status(
                "finite Schur repair input",
                schur_repair_ok,
                "The finite augmented Schur repair consequence theorem is closed.",
            ),
            "finiteRepairOperatorInputStatus": status(
                "finite repair operator input",
                repair_operator_ok,
                "The finite repair operator D_aug,N is nonnegative.",
            ),
            "finiteRepairedCoreFormInputStatus": status(
                "finite repaired core form input",
                repaired_form_ok,
                "The finite repaired Schur core form is nonnegative.",
            ),
            "fixedCoreSchurPositiveConsequenceStatus": status(
                "fixed-core Schur positivity consequence",
                ok,
                "Therefore q_N^Schur>=0 on every fixed augmented trace core.",
            ),
        },
        "finiteSchurCorePositiveClosed": ok,
        "finiteRepairedCoreFormPositiveClosed": ok,
        "fixedCoreOnly": True,
        "proof": [
            "Import the finite augmented Schur repair consequence theorem.",
            "Read off finite repaired form positivity and D_aug,N>=0.",
            "Expose only q_N^Schur>=0 for downstream fixed-core cone summation.",
        ],
        "remainingAnalyticGap": None if ok else "Finite Schur repaired core positivity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented fixed-core Schur positivity consequence theorem")
    print(f"  Schur repair input: {schur_repair_ok}")
    print(f"  D nonnegative input: {repair_operator_ok}")
    print(f"  repaired core form input: {repaired_form_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
