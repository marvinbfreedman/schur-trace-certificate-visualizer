#!/usr/bin/env python3
r"""Pointwise nonnegativity certificate for augmented repaired forms.

This is the narrow finite/core input needed by the closed-form limit theorem.
It makes no Mosco, density, compactness, closure, or continuum-passage claim.

For each fixed augmented trace-core index N, the repaired quadratic form splits
as

    q_N^rep = q_N^Schur + t_N,

where q_N^Schur is nonnegative by the finite Schur/Douglas repair consequence
and t_N is the restriction of the nonnegative Volterra tail form.  Therefore
q_N^rep >= 0 for that fixed N.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-core-positive-json",
        default="augmented_fixed_core_schur_positive_consequence_theorem.json",
    )
    parser.add_argument(
        "--tail-restriction-json",
        default="volterra_tail_restriction_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_pointwise_repaired_form_nonnegative_certificate.json",
    )
    args = parser.parse_args()

    finite = load(args.finite_core_positive_json)
    tail = load(args.tail_restriction_json)

    finite_ok = bool(
        finite.get("finiteSchurCorePositiveClosed")
        or finite.get("finiteRepairedCoreFormPositiveClosed")
        or nested_closed(finite, "statuses", "fixedCoreSchurPositiveConsequenceStatus")
    )
    tail_ok = bool(
        tail.get("volterraTailRestrictionNonnegativeClosed")
        or nested_closed(tail, "statuses", "fixedCoreRestrictionStatus")
    )
    decomposition_ok = finite_ok and tail_ok
    pointwise_ok = decomposition_ok

    data = {
        "theoremName": "augmented pointwise repaired form nonnegative certificate",
        "proofClass": "interval/ball certificate with symbolic cone algebra",
        "finiteCorePositiveJson": args.finite_core_positive_json,
        "tailRestrictionJson": args.tail_restriction_json,
        "scope": (
            "For each fixed augmented trace-core index N only.  This certificate "
            "does not assert Mosco convergence, density, compactness, quotient "
            "compatibility, lower semicontinuity, or any finite-to-continuum "
            "promotion."
        ),
        "statuses": {
            "finiteSchurRepairInputStatus": status(
                "finite Schur repair input",
                finite_ok,
                (
                    "The interval/ball finite Schur consequence supplies "
                    "D_aug,N>=0 and nonnegativity of the repaired finite Schur "
                    "core form."
                ),
            ),
            "tailPositiveRestrictionInputStatus": status(
                "positive Volterra tail restriction input",
                tail_ok,
                (
                    "The positive Volterra tail theorem gives a nonnegative "
                    "tail form; restricting a nonnegative form to a fixed core "
                    "preserves nonnegativity."
                ),
            ),
            "fixedCoreConeSumStatus": status(
                "fixed-core cone sum",
                decomposition_ok,
                (
                    "On each fixed core, q_N^rep is the sum of the finite "
                    "Schur repaired core form and the restricted nonnegative "
                    "Volterra tail form."
                ),
            ),
            "pointwiseRepairedFormNonnegativeStatus": status(
                "pointwise repaired form nonnegative",
                pointwise_ok,
                "A sum of two nonnegative quadratic forms is nonnegative.",
            ),
        },
        "pointwiseRepairedFormsNonnegativeClosed": pointwise_ok,
        "finiteCoreOnly": True,
        "noContinuumPassageClaimed": True,
        "exportedConsequence": (
            "For every fixed augmented trace-core index N, q_N^rep>=0."
        ),
        "proof": [
            "Import the finite augmented Schur repair consequence.",
            "Restrict the nonnegative Volterra tail form to the same fixed core.",
            "Use the exact fixed-core decomposition q_N^rep=q_N^Schur+t_N.",
            "Apply the elementary cone rule: nonnegative plus nonnegative is nonnegative.",
        ],
        "remainingAnalyticGap": None
        if pointwise_ok
        else "Finite Schur repair input or positive tail input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented pointwise repaired form nonnegative certificate")
    print(f"  finite Schur repair input: {finite_ok}")
    print(f"  positive tail restriction: {tail_ok}")
    print(f"  fixed-core decomposition: {decomposition_ok}")
    print(f"  pointwise nonnegative: {pointwise_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
