#!/usr/bin/env python3
r"""Terminal consequence for high-block spectral projection stability.

The high-block tail-passage theorem only needs the boolean contract

    highBlockSpectralProjectionClosed = true

or its historical alias.  The detailed spectral input theorem imports compact
compression, source gap, and abstract Riesz projection machinery; this wrapper
keeps those internals below the interface.
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
        "--spectral-input-json",
        default="high_block_spectral_projection_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_spectral_projection_consequence_theorem.json",
    )
    args = parser.parse_args()

    spectral = load(args.spectral_input_json)
    spectral_ok = bool(
        spectral.get("highBlockSpectralProjectionClosed")
        or spectral.get("compactSourceRieszProjectionConvergenceClosed")
        or nested_closed(spectral, "highBlockSpectralProjectionStatus")
    )

    data: dict[str, Any] = {
        "theoremName": "high-block spectral projection terminal consequence theorem",
        "proofClass": "symbolic identity",
        "spectralInputJson": args.spectral_input_json,
        "highBlockSpectralProjectionClosed": spectral_ok,
        "compactSourceRieszProjectionConvergenceClosed": spectral_ok,
        "sourceActiveInactiveSplittingStabilityClosed": spectral_ok,
        "statuses": {
            "highBlockSpectralProjectionConsequenceStatus": status(
                "high-block spectral projection stability consequence",
                spectral_ok,
                (
                    "The detailed high-block spectral projection input theorem "
                    "closes the compact-source Riesz projection stability "
                    "contract used by high-block tail passage."
                ),
            )
        },
        "proof": [
            "Import the detailed high-block spectral projection input theorem.",
            "Read off the closed Riesz/spectral projection stability contract.",
            "Export only the terminal high-block spectral projection consequence.",
        ],
        "remainingAnalyticGap": None
        if spectral_ok
        else "High-block spectral projection input theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block spectral projection consequence theorem")
    print(f"  spectral projection closed: {spectral_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
