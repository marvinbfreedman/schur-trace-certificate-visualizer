#!/usr/bin/env python3
r"""Narrow consequence for high-block Mosco projection convergence.

The high-block compact-compression input only needs the terminal projection
convergence flags:

    highBlockMoscoProjectionClosed = true,
    strongProjectionConvergenceClosed = true.

This file reads the full high-block Mosco projection input theorem and exposes
only those consequences.  Trace-frame recovery and the abstract projection
theorem remain below the full input theorem.
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
        "--projection-json",
        default="high_block_mosco_projection_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_mosco_projection_consequence_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.projection_json)
    mosco_ok = bool(projection.get("highBlockMoscoProjectionClosed"))
    strong_ok = bool(projection.get("strongProjectionConvergenceClosed"))
    ok = mosco_ok and strong_ok

    data = {
        "theoremName": "high-block Mosco projection consequence theorem",
        "proofClass": "symbolic identity",
        "source": "high-block Mosco projection input theorem",
        "statuses": {
            "highBlockMoscoProjectionInputStatus": status(
                "high-block Mosco projection input",
                mosco_ok,
                "The full input theorem closes the high-block Mosco projection passage.",
            ),
            "strongProjectionConvergenceInputStatus": status(
                "strong projection convergence input",
                strong_ok,
                "The full input theorem records strong convergence of the high-block projections.",
            ),
            "projectionConsequenceStatus": status(
                "projection convergence consequence",
                ok,
                "Only the terminal projection-convergence flags are exported upstream.",
            ),
        },
        "highBlockMoscoProjectionConsequenceClosed": ok,
        "highBlockMoscoProjectionClosed": mosco_ok,
        "strongProjectionConvergenceClosed": strong_ok,
        "proof": [
            "Import the full high-block Mosco projection input theorem.",
            "Export only the projection-convergence consequences used by compact compression.",
        ],
        "notExportedHere": [
            "trace-frame recovery theorem",
            "bounded sampled-trace correction construction",
            "abstract Mosco projection theorem",
        ],
        "remainingAnalyticGap": None if ok else "High-block Mosco projection input theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block Mosco projection consequence theorem")
    print(f"  projection consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
