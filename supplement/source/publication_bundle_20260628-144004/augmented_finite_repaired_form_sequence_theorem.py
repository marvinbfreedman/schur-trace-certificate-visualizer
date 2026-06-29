#!/usr/bin/env python3
r"""Finite nonnegative repaired-form sequence theorem."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-consequence-json",
        default="xi_augmented_finite_schur_repair_consequence_theorem.json",
    )
    parser.add_argument(
        "--tail-positive-json",
        default="volterra_tail_positive_form_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_finite_repaired_form_sequence_theorem.json")
    args = parser.parse_args()

    finite = load(args.finite_consequence_json)
    tail = load(args.tail_positive_json)

    finite_ok = bool(finite.get("finiteAugmentedSchurRepairClosed"))
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )
    sequence_ok = finite_ok and tail_ok

    data = {
        "theoremName": "augmented finite repaired form sequence theorem",
        "proofClass": "analytic proof",
        "finiteConsequenceJson": args.finite_consequence_json,
        "tailPositiveJson": args.tail_positive_json,
        "statement": (
            "The augmented Schur repair gives a directed sequence of "
            "nonnegative repaired quadratic forms on the finite trace-image "
            "core spaces, with the positive Volterra tail included."
        ),
        "statuses": {
            "finiteSchurConsequenceInputStatus": status(
                "finite Schur repair consequence input",
                finite_ok,
                "The finite consequence theorem supplies D_aug,N>=0 and repaired-form nonnegativity.",
            ),
            "positiveTailInputStatus": status(
                "positive Volterra tail input",
                tail_ok,
                "The tail theorem supplies the nonnegative tail part of the repaired form.",
            ),
            "finiteRepairedFormSequenceStatus": status(
                "finite repaired-form sequence",
                sequence_ok,
                (
                    "Combining the finite repair consequence with the positive "
                    "tail gives nonnegative repaired forms on every finite core space."
                ),
            ),
        },
        "finiteRepairedFormSequenceClosed": sequence_ok,
        "finiteRepairedFormsNonnegativeClosed": sequence_ok,
        "proof": [
            "Import the finite augmented Schur repair consequence.",
            "Add the separately nonnegative Volterra tail form.",
            "The sum of nonnegative finite forms remains nonnegative on each finite core space.",
        ],
        "remainingAnalyticGap": None if sequence_ok else "Finite repair consequence or tail positivity is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented finite repaired form sequence theorem")
    print(f"  finite consequence: {finite_ok}")
    print(f"  positive tail: {tail_ok}")
    print(f"  finite repaired sequence: {sequence_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
