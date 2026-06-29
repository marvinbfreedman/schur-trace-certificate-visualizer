#!/usr/bin/env python3
r"""Narrow high-block min-max tail consequence theorem.

This wrapper exposes only the final high-block min-max passage consequence
needed by the high-block tail-estimate continuum-passage theorem.
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
        "--minmax-tail-json",
        default="high_block_minmax_tail_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_minmax_tail_consequence_theorem.json",
    )
    args = parser.parse_args()

    minmax = load(args.minmax_tail_json)
    minmax_ok = bool(minmax.get("highBlockMinmaxTailPassageClosed"))
    abstract_ok = bool(minmax.get("abstractMinmaxTailPassageClosed"))
    tail_passage_ok = bool(minmax.get("tailEstimatePassesToContinuum"))
    ok = minmax_ok and abstract_ok and tail_passage_ok

    data = {
        "theoremName": "high-block min-max tail consequence theorem",
        "proofClass": "symbolic identity",
        "minmaxTailJson": args.minmax_tail_json,
        "normalizedEpsilonDelta": minmax.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": minmax.get("finiteLowMidSchurBudget"),
        "absorptionSlack": minmax.get("absorptionSlack"),
        "statement": (
            "The high-block min-max tail input theorem supplies the final "
            "inactive-tail passage consequence needed by the continuum "
            "high-block tail estimate theorem."
        ),
        "statuses": {
            "highBlockMinmaxTailPassageStatus": status(
                "high-block min-max tail passage",
                minmax_ok,
                "Imported from the high-block min-max tail input theorem.",
            ),
            "abstractMinmaxTailPassageStatus": status(
                "abstract min-max tail passage",
                abstract_ok,
                "Imported from the high-block min-max tail input theorem.",
            ),
            "tailEstimateContinuumPassageStatus": status(
                "tail estimate continuum passage",
                tail_passage_ok,
                "Imported from the high-block min-max tail input theorem.",
            ),
            "highBlockMinmaxTailConsequenceStatus": status(
                "high-block min-max tail consequence",
                ok,
                "Together these flags give the min-max tail input used upstream.",
            ),
        },
        "highBlockMinmaxTailPassageClosed": minmax_ok,
        "abstractMinmaxTailPassageClosed": abstract_ok,
        "tailEstimatePassesToContinuum": tail_passage_ok,
        "highBlockMinmaxTailConsequenceClosed": ok,
        "proof": [
            "Import the high-block min-max tail input theorem.",
            "Extract only the final min-max/tail-passage consequences used upstream.",
            "Expose these consequences to the high-block tail-estimate continuum-passage theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "High-block min-max tail passage, abstract passage, or tail continuum passage is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block min-max tail consequence theorem")
    print(f"  high-block min-max tail passage: {minmax_ok}")
    print(f"  abstract min-max tail passage: {abstract_ok}")
    print(f"  tail estimate passes to continuum: {tail_passage_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
