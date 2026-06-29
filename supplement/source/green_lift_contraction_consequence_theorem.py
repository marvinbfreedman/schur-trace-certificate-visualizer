#!/usr/bin/env python3
r"""Narrow completed Green-lift contraction consequence.

This wrapper exposes only the consequence needed by the source-range
Hardy/Green theorem:

    ||C K E|| <= 1

on the completed Volterra trace-fiber image, with the integration-by-parts
identity already passed to the closed form domain.
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
        "--green-closure-json",
        default="continuum_green_lift_closure_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="green_lift_contraction_consequence_theorem.json",
    )
    args = parser.parse_args()

    green = load(args.green_closure_json)
    closure_ok = bool(
        green.get("closedOnCompletedVolterraDomain")
        or green.get("continuumGreenLiftClosureClosed")
        or nested_closed(green, "statuses", "greenLiftContractionStatus")
    )
    ibp_ok = nested_closed(green, "statuses", "integrationByPartsClosureStatus")
    closed_kernel_ok = nested_closed(green, "statuses", "closedTraceKernelStatus")
    ok = closure_ok and ibp_ok and closed_kernel_ok

    data = {
        "theoremName": "Green-lift contraction consequence theorem",
        "proofClass": "symbolic identity",
        "greenClosureJson": args.green_closure_json,
        "statement": (
            "The continuum Green-lift closure theorem implies the compressed "
            "Hardy/Volterra Green-lift contraction on the completed trace-fiber "
            "domain."
        ),
        "statuses": {
            "greenClosureInputStatus": status(
                "completed Green-lift closure input",
                closure_ok,
                (
                    "The continuum Green-lift closure theorem proves the "
                    "compressed contraction on the completed Volterra trace "
                    "image."
                ),
            ),
            "integrationByPartsClosureInputStatus": status(
                "integration-by-parts closure input",
                ibp_ok,
                (
                    "The same theorem passes the core integration-by-parts "
                    "identity and boundary cancellation to the completed "
                    "form domain."
                ),
            ),
            "closedTraceKernelInputStatus": status(
                "closed trace-kernel input",
                closed_kernel_ok,
                "The trace kernel is closed in the completed form domain.",
            ),
            "greenLiftContractionConsequenceStatus": status(
                "Green-lift contraction consequence",
                ok,
                "Therefore ||C K E|| <= 1 on the completed trace-fiber image.",
            ),
        },
        "greenLiftContractionClosed": ok,
        "closedTraceFiberContractionClosed": ok,
        "completedDomainIntegrationByPartsClosed": ok,
        "proof": [
            "Import the continuum Green-lift closure theorem.",
            "Read off the closed trace-kernel and integration-by-parts closure statuses.",
            "Read off the compressed Green-lift contraction status.",
            "Expose only this contraction consequence to downstream source-range estimates.",
        ],
        "remainingAnalyticGap": None if ok else "Continuum Green-lift closure input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Green-lift contraction consequence theorem")
    print(f"  closure input: {closure_ok}")
    print(f"  integration-by-parts closure: {ibp_ok}")
    print(f"  closed trace kernel: {closed_kernel_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
