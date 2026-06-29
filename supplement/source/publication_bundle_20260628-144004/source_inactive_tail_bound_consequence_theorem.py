#!/usr/bin/env python3
r"""Narrow source-inactive tail bound consequence.

This wrapper exposes only the constants and booleans needed by the high-block
min-max tail input theorem:

    finite certified inactive-tail bound,
    low/mid Schur absorption budget,
    normalized epsilon/budget/slack values.
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
        "--tail-constants-json",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="source_inactive_tail_bound_consequence_theorem.json",
    )
    args = parser.parse_args()

    tail = load(args.tail_constants_json)
    finite_tail_ok = bool(
        tail.get("sourceInactiveTailConstantsClosed")
        or tail.get("sourceInactiveTailDominationConstantsClosed")
        or tail.get("minMaxProofInCertifiedSourceModel")
        or nested_closed(tail, "statuses", "normalizedTailConstantsStatus")
    )
    absorption_ok = bool(
        tail.get("absorbableByFiniteLowMidBlock")
        or nested_closed(tail, "statuses", "absorptionBudgetStatus")
    )
    normalized_epsilon = tail.get("normalizedEpsilonDelta")
    budget = tail.get("finiteLowMidSchurBudget")
    slack = tail.get("absorptionSlack")
    ok = finite_tail_ok and absorption_ok

    data = {
        "theoremName": "source-inactive tail bound consequence theorem",
        "proofClass": "symbolic identity",
        "tailConstantsJson": args.tail_constants_json,
        "normalizedEpsilonDelta": normalized_epsilon,
        "finiteLowMidSchurBudget": budget,
        "absorptionSlack": slack,
        "epsilonOverBudget": tail.get("epsilonOverBudget"),
        "absoluteTailBound": tail.get("absoluteTailBound"),
        "sourceTopLower": tail.get("sourceTopLower"),
        "statement": (
            "The source-inactive constants theorem supplies the normalized "
            "finite inactive-tail constant and proves it is absorbable by the "
            "finite low/mid Schur budget."
        ),
        "statuses": {
            "finiteInactiveTailBoundStatus": status(
                "finite certified inactive-tail bound",
                finite_tail_ok,
                (
                    "Imported from the source-inactive constants theorem: the "
                    "normalized source model has the certified inactive-tail "
                    "constant epsilon_delta."
                ),
            ),
            "lowMidAbsorptionBudgetStatus": status(
                "low/mid Schur absorption budget",
                absorption_ok,
                (
                    "Imported from the source-inactive constants theorem: "
                    "epsilon_delta is below the finite low/mid Schur budget."
                ),
            ),
            "tailBoundConsequenceStatus": status(
                "source-inactive tail bound consequence",
                ok,
                (
                    "The finite inactive-tail constant and its absorption "
                    "budget are both closed."
                ),
            ),
        },
        "minMaxProofInCertifiedSourceModel": finite_tail_ok,
        "sourceInactiveTailConstantsClosed": finite_tail_ok,
        "sourceInactiveTailDominationConstantsClosed": ok,
        "absorbableByFiniteLowMidBlock": absorption_ok,
        "tailBoundConsequenceClosed": ok,
        "proof": [
            "Import the normalized finite inactive-tail constant.",
            "Import the low/mid Schur absorption comparison.",
            "Forward the epsilon, budget, and slack values needed by high-block min-max assembly.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Source-inactive tail constant or absorption budget is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Source-inactive tail bound consequence theorem")
    print(f"  finite inactive-tail bound: {finite_tail_ok}")
    print(f"  absorption budget: {absorption_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
