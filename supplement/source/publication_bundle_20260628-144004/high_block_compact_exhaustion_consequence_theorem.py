#!/usr/bin/env python3
r"""Symbolic consequence wrapper for high-block compact exhaustion."""

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
        "--tail-passage-json",
        default="high_block_tail_estimate_continuum_passage_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_compact_exhaustion_consequence_theorem.json",
    )
    args = parser.parse_args()

    passage = load(args.tail_passage_json)
    passage_ok = bool(passage.get("tailEstimatePassesToContinuum"))

    data = {
        "theoremName": "high-block compact exhaustion consequence identity",
        "proofClass": "symbolic identity",
        "normalizedEpsilonDelta": passage.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": passage.get("finiteLowMidSchurBudget"),
        "absorptionSlack": passage.get("absorptionSlack"),
        "statuses": {
            "tailPassageInputStatus": status(
                "high-block tail passage input",
                passage_ok,
                "The high-block tail-passage theorem proves the continuum estimate.",
            ),
            "compactExhaustionConsequenceStatus": status(
                "compact exhaustion consequence",
                passage_ok,
                (
                    "The statement 'high-block compact exhaustion is closed' "
                    "is exactly the continuum tail-passage theorem in the "
                    "high-block source model."
                ),
            ),
        },
        "conditionalHighBlockExhaustionClosed": passage_ok,
        "tailEstimatePassesToContinuum": passage_ok,
        "proof": [
            "Import the high-block tail estimate continuum passage theorem.",
            "Expose its conclusion under the historical compact-exhaustion names.",
        ],
        "remainingAnalyticGap": None if passage_ok else "High-block tail passage input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block compact exhaustion consequence theorem")
    print(f"  tail passage input: {passage_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
