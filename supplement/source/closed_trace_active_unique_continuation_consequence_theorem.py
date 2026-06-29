#!/usr/bin/env python3
r"""Terminal consequence of closed-trace active unique continuation."""

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
        "--unique-continuation-json",
        default="closed_trace_active_unique_continuation_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="closed_trace_active_unique_continuation_consequence_theorem.json",
    )
    args = parser.parse_args()

    theorem = load(args.unique_continuation_json)
    closed = bool(theorem.get("closedTraceActiveUniqueContinuationClosed"))

    data = {
        "theoremName": "closed-trace active unique-continuation consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "closedTraceActiveUniqueContinuationConsequenceStatus": status(
                "closed-trace active unique-continuation consequence",
                closed,
                (
                    "The closed-trace active unique-continuation theorem "
                    "proves Lambda_a(f)=0 on the trace interval implies "
                    "E_active f=0 in the completed closed trace domain."
                ),
            )
        },
        "closedTraceActiveUniqueContinuationClosed": closed,
        "activeSourceRowsVanishOnClosedTraceKernel": closed,
        "remainingAnalyticGap": None
        if closed
        else "Close closed_trace_active_unique_continuation_theorem.py.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed-trace active unique-continuation consequence theorem")
    print(f"  unique continuation closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
