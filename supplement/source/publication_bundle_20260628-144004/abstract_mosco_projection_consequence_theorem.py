#!/usr/bin/env python3
r"""Narrow abstract Mosco projection convergence consequence.

This wrapper exposes only the conclusion needed by model-specific high-block
code:

    Mosco convergence of closed subspaces gives strong convergence of their
    Hilbert orthogonal projections.
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
        "--abstract-projection-json",
        default="abstract_mosco_projection_convergence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_mosco_projection_consequence_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.abstract_projection_json)
    strong_ok = bool(
        projection.get("strongProjectionConvergenceClosed")
        or projection.get("moscoProjectionConvergenceClosed")
        or nested_closed(projection, "statuses", "projectionStrongConvergenceStatus")
    )

    data = {
        "theoremName": "abstract Mosco projection consequence theorem",
        "proofClass": "symbolic identity",
        "abstractProjectionJson": args.abstract_projection_json,
        "statement": (
            "The abstract Mosco projection theorem implies strong convergence "
            "of the Hilbert orthogonal projections for Mosco-convergent closed "
            "subspaces."
        ),
        "statuses": {
            "abstractMoscoProjectionInputStatus": status(
                "abstract Mosco projection theorem input",
                strong_ok,
                "The abstract theorem proves strong projection convergence.",
            ),
            "strongProjectionConvergenceStatus": status(
                "strong projection convergence consequence",
                strong_ok,
                (
                    "This is the only abstract functional-analytic conclusion "
                    "needed by the high-block Mosco projection input theorem."
                ),
            ),
        },
        "strongProjectionConvergenceClosed": strong_ok,
        "moscoProjectionConvergenceClosed": strong_ok,
        "proof": [
            "Import the abstract Mosco projection convergence theorem.",
            "Read off strong Hilbert projection convergence.",
            "Expose only that consequence to model-specific high-block code.",
        ],
        "remainingAnalyticGap": None
        if strong_ok
        else "Abstract Mosco projection convergence theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract Mosco projection consequence theorem")
    print(f"  strong projection convergence: {strong_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
