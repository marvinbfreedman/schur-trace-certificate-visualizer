#!/usr/bin/env python3
r"""Narrow active-trace control consequence for source-range estimates.

This wrapper exposes only the consequence needed by the hardened
source-range Hardy/Green theorem:

    f in H_M cap ker R_global  =>  E_active f = 0.

The detailed proof remains in the synchronized active-range theorem and the
closed-trace active unique-continuation theorem.
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


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--active-range-json",
        default="synchronized_active_range_interval_consequence_theorem.json",
    )
    parser.add_argument(
        "--unique-continuation-json",
        default="closed_trace_active_unique_continuation_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="active_trace_control_consequence_theorem.json",
    )
    args = parser.parse_args()

    active = load(args.active_range_json)
    unique = load(args.unique_continuation_json)

    active_range_closed = closed(active, "activeRangeInclusionStatus")
    active_annihilation_closed = closed(active, "closedTraceActiveAnnihilationStatus")
    unique_closed = bool(unique.get("closedTraceActiveUniqueContinuationClosed"))
    ok = active_range_closed and active_annihilation_closed and unique_closed

    data = {
        "theoremName": "active trace control consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "On the completed high closed-trace block, the active source rows "
            "factor through the global trace map and hence vanish on ker R_global."
        ),
        "statuses": {
            "synchronizedActiveRangeInputStatus": status(
                "synchronized active range input",
                active_range_closed and active_annihilation_closed,
                (
                    "The synchronized active-range theorem supplies active "
                    "range inclusion and closed-trace active annihilation."
                ),
            ),
            "closedTraceUniqueContinuationInputStatus": status(
                "closed-trace active unique-continuation input",
                unique_closed,
                (
                    "The closed-trace active unique-continuation theorem "
                    "proves Lambda_a(f)=0 on I implies E_active f=0."
                ),
            ),
            "activeTraceControlConsequenceStatus": status(
                "active trace control consequence",
                ok,
                (
                    "Combining the active-range inclusion with closed-trace "
                    "unique continuation gives R_global f=0 => E_active f=0."
                ),
            ),
        },
        "activeTraceControlClosed": ok,
        "activeSourceRowsVanishOnClosedTraceKernel": ok,
        "proof": [
            "Use the closed-trace active unique-continuation theorem on the interval trace family.",
            "Use the synchronized active-range theorem to identify the active source plane with the closed trace range.",
            "Conclude that every active source row annihilates H_M cap ker R_global.",
        ],
        "remainingAnalyticGap": None if ok else "Active range or unique-continuation input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Active trace control consequence theorem")
    print(f"  synchronized active range input: {active_range_closed and active_annihilation_closed}")
    print(f"  unique-continuation input: {unique_closed}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
