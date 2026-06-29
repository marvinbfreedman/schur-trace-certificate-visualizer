#!/usr/bin/env python3
r"""Terminal consequence for the high-block exhaustion theorem.

The quotient Schur assembly only needs the consequence

    tailEstimatePassesToContinuum = true

plus the scalar tail/budget constants.  It does not need to import the full
high-block exhaustion ledger and thereby pull the detailed Mosco, compactness,
trace-frame, and source-operator proof chain into every upstream audit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compact-consequence-json",
        default="high_block_compact_exhaustion_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_exhaustion_consequence_theorem.json",
    )
    args = parser.parse_args()

    compact = load(args.compact_consequence_json)
    tail_ok = bool(compact.get("tailEstimatePassesToContinuum"))
    conditional_ok = bool(compact.get("conditionalHighBlockExhaustionClosed"))
    theorem_closed = tail_ok and conditional_ok

    data: dict[str, Any] = {
        "theoremName": "high-block exhaustion terminal consequence theorem",
        "proofClass": "symbolic identity",
        "normalizedEpsilonDelta": compact.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": compact.get("finiteLowMidSchurBudget"),
        "absorptionSlack": compact.get("absorptionSlack"),
        "tailEstimatePassesToContinuum": theorem_closed,
        "conditionalHighBlockExhaustionClosed": theorem_closed,
        "statuses": {
            "compactExhaustionInputStatus": status(
                "compact-exhaustion consequence input",
                conditional_ok,
                (
                    "The compact-exhaustion consequence supplies the closed "
                    "high-block exhaustion implication."
                ),
            ),
            "tailPassageConsequenceStatus": status(
                "continuum tail-passage terminal consequence",
                theorem_closed,
                (
                    "The continuum high-block tail estimate may be used by "
                    "quotient Schur assembly without importing the full "
                    "high-block proof ledger."
                ),
            ),
        },
        "proof": [
            "Import the compact-exhaustion consequence theorem.",
            "Read off the terminal continuum tail-passage conclusion.",
            "Export only the high-block exhaustion contract used by quotient assembly.",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "High-block compact-exhaustion consequence is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block exhaustion terminal consequence theorem")
    print(f"  compact consequence input: {conditional_ok}")
    print(f"  tail estimate passes to continuum: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
