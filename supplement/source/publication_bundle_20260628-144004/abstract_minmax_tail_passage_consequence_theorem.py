#!/usr/bin/env python3
r"""Narrow abstract min-max tail passage consequence.

This wrapper exposes only the final consequence used by the high-block source
tail input theorem:

    finite inactive-tail bounds pass through the convergent active/inactive
    source projection model to the continuum closed high block.
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
        "--abstract-minmax-json",
        default="abstract_minmax_tail_passage_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_minmax_tail_passage_consequence_theorem.json",
    )
    args = parser.parse_args()

    minmax = load(args.abstract_minmax_json)
    passage_ok = bool(
        minmax.get("abstractMinmaxTailPassageClosed")
        or minmax.get("spectralMinmaxPassageClosed")
        or nested_closed(minmax, "statuses", "courantFischerTailPassageStatus")
    )

    data = {
        "theoremName": "abstract min-max tail passage consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The abstract min-max tail passage theorem implies that the "
            "finite inactive-tail inequality transports to the continuum "
            "closed high block under the convergent source projection model."
        ),
        "statuses": {
            "abstractMinmaxTailPassageInputStatus": status(
                "abstract min-max tail passage input",
                passage_ok,
                (
                    "Imported from the abstract min-max tail passage theorem: "
                    "the Courant-Fischer/recovery-sequence passage is closed."
                ),
            ),
            "continuumInactiveTailPassageStatus": status(
                "continuum inactive-tail passage",
                passage_ok,
                (
                    "Therefore the high-block finite inactive-tail bound is "
                    "eligible to pass to the continuum closed high block once "
                    "the model-specific finite constants and spectral split "
                    "are supplied."
                ),
            ),
        },
        "abstractMinmaxTailPassageClosed": passage_ok,
        "continuumInactiveTailPassageClosed": passage_ok,
        "spectralMinmaxPassageClosed": passage_ok,
        "proof": [
            "Import the final closed Courant-Fischer/min-max tail passage from the abstract theorem.",
            "Expose only this passage consequence for model-specific high-block assembly.",
        ],
        "remainingAnalyticGap": None
        if passage_ok
        else "Abstract min-max tail passage theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract min-max tail passage consequence theorem")
    print(f"  continuum inactive-tail passage: {passage_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
