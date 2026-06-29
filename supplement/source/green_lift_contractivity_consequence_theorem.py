#!/usr/bin/env python3
r"""Slim consequence for Green-lift contractivity.

The continuum Green-lift closure theorem only consumes the closed contraction
statement

    ||C K E|| <= 1

on the completed Green-minimizer trace image.  This wrapper exports that
single consequence from the full closed-form contractivity theorem, leaving the
Euler boundary/minimality and Hardy multiplier ingredients below it in the
audit graph.
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
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contractivity-json",
        default="green_lift_contractivity_form_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="green_lift_contractivity_consequence_theorem.json",
    )
    args = parser.parse_args()

    contractivity = load(args.contractivity_json)
    contraction_ok = bool(
        contractivity.get("greenLiftContractionClosed")
        or contractivity.get("closedTraceFiberContractionClosed")
        or nested_closed(contractivity, "statuses", "compressedGreenLiftContractionStatus")
    )

    data = {
        "theoremName": "Green-lift contractivity consequence theorem",
        "proofClass": "symbolic consequence",
        "contractivityJson": args.contractivity_json,
        "statuses": {
            "compressedGreenLiftContractionStatus": status(
                "compressed Green-lift contraction",
                contraction_ok,
                (
                    "The imported closed-form Green-lift contractivity theorem "
                    "proves ||C K E||<=1 on the completed Green-minimizer "
                    "trace image."
                ),
            )
        },
        "closedTraceFiberContractionClosed": contraction_ok,
        "greenLiftContractionClosed": contraction_ok,
        "proof": [
            "Import the closed-form Green-lift contractivity theorem.",
            "Expose only the completed trace-fiber contraction needed by the closure theorem.",
            "Keep boundary minimality and Hardy multiplier inputs below the full contractivity theorem.",
        ],
        "remainingAnalyticGap": None
        if contraction_ok
        else "Full Green-lift contractivity theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Green-lift contractivity consequence theorem")
    print(f"  contraction closed: {contraction_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
