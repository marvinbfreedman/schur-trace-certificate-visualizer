#!/usr/bin/env python3
r"""Narrow consequence for high-block source-operator compactness.

The high-block compact-compression input only needs the terminal compactness
statement:

    compactSourceOperatorClosed = true.

This file reads the full high-block source operator compactness theorem and
exports only that consequence.  The Hardy/Green estimate, representer
continuity, and compact Bochner-integral proof remain below the full theorem.
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
        "--source-compactness-json",
        default="high_block_source_operator_compactness_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_source_operator_compactness_consequence_theorem.json",
    )
    args = parser.parse_args()

    source = load(args.source_compactness_json)
    compact_ok = bool(source.get("compactSourceOperatorClosed"))
    high_block_ok = bool(source.get("highBlockSourceOperatorCompactnessClosed"))
    ok = compact_ok and high_block_ok

    data = {
        "theoremName": "high-block source operator compactness consequence theorem",
        "proofClass": "symbolic identity",
        "source": "high-block source operator compactness theorem",
        "statuses": {
            "compactSourceOperatorInputStatus": status(
                "compact source operator input",
                compact_ok,
                "The full source compactness theorem closes compactness of S=E^*E.",
            ),
            "highBlockSourceCompactnessInputStatus": status(
                "high-block source compactness input",
                high_block_ok,
                "The full theorem records compactness on the completed closed-trace high block.",
            ),
            "sourceCompactnessConsequenceStatus": status(
                "source compactness consequence",
                ok,
                "Only the terminal compact-source flag is exported upstream.",
            ),
        },
        "highBlockSourceOperatorCompactnessConsequenceClosed": ok,
        "highBlockSourceOperatorCompactnessClosed": high_block_ok,
        "compactSourceOperatorClosed": compact_ok,
        "proof": [
            "Import the full high-block source operator compactness theorem.",
            "Export only compactness of the high-block source operator used by compact compression.",
        ],
        "notExportedHere": [
            "closed-trace Hardy/Green estimate",
            "uniform Riesz representer construction",
            "finite-rank source approximation",
        ],
        "remainingAnalyticGap": None if ok else "High-block source operator compactness theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block source operator compactness consequence theorem")
    print(f"  source compactness consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
