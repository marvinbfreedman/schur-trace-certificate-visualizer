#!/usr/bin/env python3
r"""Narrow augmented trace limit-identification consequence.

This wrapper exposes only the consequence needed by the augmented trace liminf
representative compactness theorem:

    closed augmented trace convergence identifies the weak graph limit's
    augmented trace coordinate.

The theta-tail, Mellin-boundary, and Green-closure inputs remain inside
``xi_augmented_trace_convergence_theorem.json``.
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
        "--trace-convergence-json",
        default="xi_augmented_trace_convergence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_trace_limit_identification_consequence_theorem.json",
    )
    args = parser.parse_args()

    trace = load(args.trace_convergence_json)
    closed_trace_ok = bool(
        trace.get("closedTraceConvergenceClosed")
        or nested_closed(trace, "statuses", "closedTraceConvergenceStatus")
    )
    graph_trace_ok = bool(trace.get("augmentedTraceGraphConvergenceClosed"))
    limit_identification_ok = closed_trace_ok and graph_trace_ok

    data = {
        "theoremName": "Xi augmented trace limit-identification consequence theorem",
        "proofClass": "symbolic identity",
        "traceConvergenceJson": args.trace_convergence_json,
        "statement": (
            "Strong graph-dual convergence of the finite augmented traces "
            "R_aug,N to the closed trace R_aug identifies the augmented trace "
            "coordinate of weak graph limits in the liminf representative "
            "compactness argument."
        ),
        "statuses": {
            "closedAugmentedTraceConvergenceInputStatus": status(
                "closed augmented trace convergence input",
                closed_trace_ok,
                (
                    "Imported from the Xi augmented trace convergence theorem: "
                    "R_aug,N converges to the closed augmented trace R_aug in "
                    "the transported graph-dual topology."
                ),
            ),
            "augmentedTraceGraphConvergenceInputStatus": status(
                "augmented trace graph convergence input",
                graph_trace_ok,
                (
                    "The same theorem exports convergence in the completed "
                    "augmented trace graph topology."
                ),
            ),
            "traceLimitIdentificationStatus": status(
                "trace limit identification",
                limit_identification_ok,
                (
                    "For bounded representative sequences with weak graph "
                    "limits, the strong graph-dual trace convergence identifies "
                    "the limiting augmented trace coordinate."
                ),
            ),
        },
        "closedTraceConvergenceClosed": closed_trace_ok,
        "augmentedTraceGraphConvergenceClosed": graph_trace_ok,
        "traceLimitIdentificationClosed": limit_identification_ok,
        "proof": [
            "Import closed augmented trace convergence.",
            "Import augmented trace graph convergence.",
            "Apply these convergences to the weak graph limit supplied by the liminf representative compactness argument.",
            "Expose only trace-limit identification to downstream quotient-liminf code.",
        ],
        "remainingAnalyticGap": None
        if limit_identification_ok
        else "Closed augmented trace convergence input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Xi augmented trace limit-identification consequence theorem")
    print(f"  closed trace convergence: {closed_trace_ok}")
    print(f"  graph trace convergence: {graph_trace_ok}")
    print(f"  trace limit identification: {limit_identification_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
