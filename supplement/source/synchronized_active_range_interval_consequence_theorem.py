#!/usr/bin/env python3
r"""Terminal consequence for synchronized active-range interval control.

The high-block tail-passage theorem only needs the active-range annihilation
contract and the source-inactive tail status.  This wrapper exports those
terminal facts without importing the full endpoint/Krawczyk active-range proof
ledger upstream.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def nested_closed(data: dict[str, Any], key: str) -> bool:
    item = data.get(key, {})
    if not item and isinstance(data.get("statuses"), dict):
        item = data["statuses"].get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


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
        "--active-range-json",
        default="synchronized_active_range_interval_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="synchronized_active_range_interval_consequence_theorem.json",
    )
    args = parser.parse_args()

    active = load(args.active_range_json)
    endpoint_ok = nested_closed(active, "endpointGreenBvpSolvabilityStatus")
    annihilation_ok = nested_closed(active, "closedTraceActiveAnnihilationStatus")
    active_range_ok = nested_closed(active, "activeRangeInclusionStatus")
    inactive_tail_ok = nested_closed(active, "sourceInactiveTailDominationStatus")
    combined_ok = nested_closed(active, "fullContinuumCombinedSourceBoundStatus")
    theorem_closed = active_range_ok and annihilation_ok and inactive_tail_ok

    data: dict[str, Any] = {
        "theoremName": "synchronized active-range interval terminal consequence theorem",
        "proofClass": "symbolic identity",
        "activeRangeJson": args.active_range_json,
        "endpointGreenBvpSolvabilityStatus": status(
            "endpoint Green BVP solvability consequence",
            endpoint_ok,
            "The full active-range interval theorem closes endpoint Green BVP solvability.",
        ),
        "closedTraceActiveAnnihilationStatus": status(
            "closed trace annihilates active source",
            annihilation_ok,
            "The full active-range interval theorem proves ker R_global annihilates E_active.",
        ),
        "activeRangeInclusionStatus": status(
            "E_active in closure Range(R_global^*)",
            active_range_ok,
            "The full active-range interval theorem proves the Hilbert annihilator/range inclusion.",
        ),
        "sourceInactiveTailDominationStatus": status(
            "source-inactive min-max tail domination",
            inactive_tail_ok,
            "The full active-range interval theorem imports the source-inactive tail consequence.",
        ),
        "fullContinuumCombinedSourceBoundStatus": status(
            "full continuum active range plus inactive tail",
            combined_ok,
            "The full active-range interval theorem combines active annihilation with the inactive tail bound.",
        ),
        "synchronizedActiveRangeIntervalConsequenceClosed": theorem_closed,
        "proof": [
            "Import the full synchronized active-range interval theorem.",
            "Read off active annihilation, active range inclusion, and inactive-tail domination.",
            "Export only the terminal active/inactive source contract used by high-block tail passage.",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "One active-range interval consequence status is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Synchronized active-range interval consequence theorem")
    print(f"  active range inclusion: {active_range_ok}")
    print(f"  closed trace annihilation: {annihilation_ok}")
    print(f"  inactive tail domination: {inactive_tail_ok}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
