#!/usr/bin/env python3
r"""Narrow high-block active/inactive source split consequence.

This wrapper exposes only the source split consequence needed by the high-block
min-max tail input theorem:

    the active two-dimensional source subspace and inactive complement are
    stable in the high-block source model.
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
        "--spectral-json",
        default="high_block_spectral_projection_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_active_source_split_consequence_theorem.json",
    )
    args = parser.parse_args()

    spectral = load(args.spectral_json)
    split_ok = bool(
        spectral.get("highBlockSpectralProjectionClosed")
        or spectral.get("sourceActiveInactiveSplittingStabilityClosed")
        or nested_closed(spectral, "statuses", "highBlockSpectralProjectionStatus")
    )

    data = {
        "theoremName": "high-block active source split consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The high-block spectral projection input theorem identifies a "
            "stable active top-two source subspace and the corresponding "
            "inactive complement for the high-block source model."
        ),
        "statuses": {
            "highBlockSpectralProjectionInputStatus": status(
                "high-block spectral projection input",
                split_ok,
                (
                    "Imported from the high-block spectral projection theorem: "
                    "compact compression, source gap, and abstract Riesz "
                    "projection stability close the high-block source split."
                ),
            ),
            "activeInactiveSourceSplitStatus": status(
                "stable active/inactive source split",
                split_ok,
                (
                    "The active top-two source subspace and inactive "
                    "orthogonal complement are stable for the high-block "
                    "min-max tail comparison."
                ),
            ),
        },
        "highBlockSpectralProjectionClosed": split_ok,
        "activeInactiveSourceSplitClosed": split_ok,
        "sourceActiveInactiveSplittingStabilityClosed": split_ok,
        "proof": [
            "Import the high-block spectral projection theorem.",
            "Read off stable active/inactive source splitting.",
            "Expose only this split consequence to the high-block min-max tail input theorem.",
        ],
        "remainingAnalyticGap": None
        if split_ok
        else "High-block spectral projection input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block active source split consequence theorem")
    print(f"  active/inactive source split: {split_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
