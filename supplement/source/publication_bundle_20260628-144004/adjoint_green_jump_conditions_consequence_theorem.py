#!/usr/bin/env python3
r"""Terminal consequence of the adjoint Green jump-conditions theorem."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {"label": label, "closed": closed, "status": "closed" if closed else "open", "reason": reason}
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jump-json", default="adjoint_green_jump_conditions_theorem.json")
    parser.add_argument("--json-out", default="adjoint_green_jump_conditions_consequence_theorem.json")
    args = parser.parse_args()

    jump = load(args.jump_json)
    closed_statements = jump.get("closedStatements", {})
    theorem_closed = bool(jump.get("theoremClosed"))
    interior_solved = bool(closed_statements.get("interiorSourceSolved"))
    closed = theorem_closed and interior_solved

    data = {
        "theoremName": "adjoint Green jump conditions consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "jumpConditionsConsequenceStatus": status(
                "adjoint Green jump conditions consequence",
                closed,
                (
                    "The distributional adjoint Green jump theorem supplies "
                    "the interior source jump law and solves the interior "
                    "source equation used by the endpoint range theorem."
                ),
                blocker=None if closed else "Close adjoint_green_jump_conditions_theorem.py.",
            )
        },
        "theoremClosed": closed,
        "closedStatements": {"interiorSourceSolved": interior_solved},
        "interiorSourceSolvedClosed": interior_solved,
        "remainingAnalyticGap": None if closed else "Close the full jump-conditions theorem.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Adjoint Green jump conditions consequence theorem")
    print(f"  interior source solved: {interior_solved}")
    print(f"  theorem closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
