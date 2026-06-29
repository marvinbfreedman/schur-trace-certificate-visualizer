#!/usr/bin/env python3
r"""Narrow consequence of the closed-LSC nonnegative cone principle."""

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
        "--cone-principle-json",
        default="closed_lsc_nonnegative_cone_principle_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="closed_lsc_nonnegative_cone_consequence_theorem.json",
    )
    args = parser.parse_args()

    cone = load(args.cone_principle_json)
    cone_ok = bool(cone.get("closedLscNonnegativeConePrincipleClosed"))

    data = {
        "theoremName": "closed LSC nonnegative cone consequence theorem",
        "proofClass": "symbolic identity",
        "conePrincipleJson": args.cone_principle_json,
        "statement": (
            "The closed-LSC nonnegative cone principle supplies the abstract "
            "cone-closure consequence needed by the augmented nonnegative "
            "closed-form limit theorem."
        ),
        "statuses": {
            "closedLscNonnegativeConePrincipleStatus": status(
                "closed LSC nonnegative cone principle",
                cone_ok,
                "Imported from the abstract closed-LSC nonnegative cone principle theorem.",
            ),
            "closedLscNonnegativeConeConsequenceStatus": status(
                "closed LSC nonnegative cone consequence",
                cone_ok,
                (
                    "A closed lower-semicontinuous envelope of nonnegative "
                    "core forms remains nonnegative."
                ),
            ),
        },
        "closedLscNonnegativeConePrincipleClosed": cone_ok,
        "closedLscNonnegativeConeConsequenceClosed": cone_ok,
        "proof": [
            "Import the abstract closed-LSC nonnegative cone principle.",
            "Extract only the cone-closure consequence used upstream.",
        ],
        "remainingAnalyticGap": None if cone_ok else "Closed-LSC cone principle input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed LSC nonnegative cone consequence theorem")
    print(f"  cone principle: {cone_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
