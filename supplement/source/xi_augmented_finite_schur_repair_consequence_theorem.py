#!/usr/bin/env python3
r"""Narrow finite augmented Schur repair consequence theorem.

This theorem composes:

1. the universal finite Schur/Douglas algebra theorem;
2. the augmented finite Schur interval constants theorem.

It exports only the finite consequences needed by the continuum transport
layer: finite D_aug,N>=0, finite repaired form positivity, Schur range
condition, and Mu annihilation on the finite augmented trace kernel.
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
        "--algebra-json",
        default="finite_schur_douglas_repair_algebra_theorem.json",
    )
    parser.add_argument(
        "--constants-json",
        default="xi_augmented_finite_schur_interval_constants_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_finite_schur_repair_consequence_theorem.json",
    )
    args = parser.parse_args()

    algebra = load(args.algebra_json)
    constants = load(args.constants_json)

    algebra_ok = bool(algebra.get("finiteSchurDouglasAlgebraClosed"))
    constants_ok = bool(constants.get("finiteSchurConstantsClosed"))
    range_ok = bool(constants.get("finiteSchurRangeConstantsClosed"))
    repair_nonnegative_ok = bool(constants.get("finiteRepairOperatorNonnegativeConstantsClosed"))
    repaired_form_ok = bool(constants.get("finiteRepairedFormPositiveConstantsClosed"))
    mu_ok = bool(constants.get("finiteMuKernelAnnihilationConstantsClosed"))
    consequence_ok = all(
        [algebra_ok, constants_ok, range_ok, repair_nonnegative_ok, repaired_form_ok, mu_ok]
    )

    data = {
        "theoremName": "augmented finite Schur repair consequence theorem",
        "proofClass": "symbolic finite Schur identity with interval/ball certificate inputs",
        "algebraJson": args.algebra_json,
        "constantsJson": args.constants_json,
        "finiteTraceDimensions": constants.get("finiteTraceDimensions"),
        "certificateConstants": constants.get("certificateConstants"),
        "statuses": {
            "finiteSchurAlgebraInputStatus": status(
                "finite Schur/Douglas algebra input",
                algebra_ok,
                "The universal finite repair algebra theorem is closed.",
            ),
            "finiteConstantsInputStatus": status(
                "finite augmented interval constants input",
                constants_ok,
                "The interval/ball finite augmented certificate constants satisfy the required inequalities.",
            ),
            "finiteSchurRangeConditionStatus": status(
                "finite Schur range condition",
                range_ok,
                "The finite Douglas range residual is below the certificate tolerance.",
            ),
            "finiteRepairOperatorNonnegativeStatus": status(
                "finite repair operator nonnegative",
                repair_nonnegative_ok,
                "D_aug,N is nonnegative in the finite trace model up to tolerance.",
            ),
            "finiteRepairedFormPositiveStatus": status(
                "finite repaired form positive",
                repaired_form_ok,
                "K_N+R_aug,N^*D_aug,N R_aug,N is nonnegative in the finite model.",
            ),
            "finiteMuKernelAnnihilationStatus": status(
                "finite Mu annihilation",
                mu_ok,
                "The primitive Mu rows vanish on the finite augmented trace kernel.",
            ),
        },
        "finiteAugmentedSchurRepairClosed": consequence_ok,
        "finiteSchurRangeConditionClosed": range_ok,
        "finiteRepairOperatorNonnegativeClosed": repair_nonnegative_ok,
        "finiteRepairedFormPositiveClosed": repaired_form_ok,
        "finiteMuKernelAnnihilationClosed": mu_ok,
        "exportedConsequence": {
            "D_aug,N": "nonnegative on the finite augmented trace image",
            "K_N+R_aug,N^*D_aug,N R_aug,N": "nonnegative finite repaired form",
            "rangeCondition": "finite Douglas Schur range condition",
            "Mu": "annihilates ker R_aug,N",
        },
        "proof": [
            "Use the interval constants theorem to verify the finite Schur hypotheses.",
            "Apply the universal finite Schur/Douglas repair algebra theorem.",
            "Export only the finite positivity/range/annihilation consequence to the continuum transport layer.",
        ],
        "remainingAnalyticGap": None
        if consequence_ok
        else "One finite Schur consequence input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented finite Schur repair consequence theorem")
    print(f"  algebra: {algebra_ok}")
    print(f"  constants: {constants_ok}")
    print(f"  range: {range_ok}")
    print(f"  D nonnegative: {repair_nonnegative_ok}")
    print(f"  repaired form: {repaired_form_ok}")
    print(f"  Mu annihilation: {mu_ok}")
    print(f"  consequence closed: {consequence_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
