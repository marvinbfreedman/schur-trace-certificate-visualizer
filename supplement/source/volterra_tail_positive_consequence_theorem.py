#!/usr/bin/env python3
r"""Narrow completed-domain consequence of Volterra tail positivity.

This wrapper exposes only the completed-domain nonnegative Volterra tail form
needed by the augmented repair transport-limit theorem.
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
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tail-positive-json",
        default="volterra_tail_positive_form_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="volterra_tail_positive_consequence_theorem.json",
    )
    args = parser.parse_args()

    tail = load(args.tail_positive_json)
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )

    data = {
        "theoremName": "completed Volterra tail positivity consequence theorem",
        "proofClass": "symbolic identity",
        "tailPositiveJson": args.tail_positive_json,
        "omegaIndependent": bool(tail.get("omegaIndependent", True)),
        "statement": (
            "The positive Volterra tail form theorem supplies the completed-domain "
            "nonnegative diagonal Volterra tail consequence used by the augmented "
            "repair transport-limit theorem."
        ),
        "statuses": {
            "volterraTailPositiveFormStatus": status(
                "positive completed Volterra tail form",
                tail_ok,
                "Imported from the positive Volterra tail form theorem.",
            ),
            "volterraTailPositiveConsequenceStatus": status(
                "completed Volterra tail positivity consequence",
                tail_ok,
                (
                    "The completed-domain diagonal Volterra tail form is "
                    "nonnegative and may be added to the transported augmented "
                    "repair form."
                ),
            ),
        },
        "volterraTailPositiveFormClosed": tail_ok,
        "completedVolterraTailPositiveConsequenceClosed": tail_ok,
        "proof": [
            "Import the positive Volterra tail form theorem.",
            "Extract only the completed-domain nonnegative tail-form consequence.",
            "Expose that consequence to the augmented repair transport-limit theorem.",
        ],
        "remainingAnalyticGap": None if tail_ok else "Volterra tail positivity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Completed Volterra tail positivity consequence theorem")
    print(f"  Volterra tail positive form: {tail_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
