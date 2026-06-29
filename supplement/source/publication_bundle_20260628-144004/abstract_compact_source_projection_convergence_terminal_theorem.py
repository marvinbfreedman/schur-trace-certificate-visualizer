#!/usr/bin/env python3
r"""Terminal projection-convergence consequence for compact source operators."""

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
        default="abstract_compact_source_spectral_projection_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_compact_source_projection_convergence_terminal_theorem.json",
    )
    args = parser.parse_args()

    spectral = load(args.spectral_json)
    riesz_ok = bool(
        spectral.get("compactSourceRieszProjectionConvergenceClosed")
        or nested_closed(spectral, "statuses", "rieszProjectionStabilityStatus")
    )
    splitting_ok = bool(
        spectral.get("sourceActiveInactiveSplittingStabilityClosed")
        or nested_closed(spectral, "statuses", "activeInactiveSplittingStabilityStatus")
    )
    closed = riesz_ok and splitting_ok

    data = {
        "theoremName": "abstract compact source projection convergence terminal theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "activeRieszProjectionConvergenceStatus": status(
                "active Riesz projection convergence terminal fact",
                closed,
                "Imported from the compact-source spectral projection theorem.",
            ),
            "activeInactiveSplittingStabilityStatus": status(
                "active/inactive splitting stability terminal fact",
                splitting_ok,
                "Imported from the compact-source spectral projection theorem.",
            ),
        },
        "activeRieszProjectionConvergenceClosed": closed,
        "inactiveProjectionConvergenceClosed": closed,
        "compactSourceRieszProjectionConvergenceClosed": closed,
        "sourceActiveInactiveSplittingStabilityClosed": splitting_ok,
        "remainingAnalyticGap": None
        if closed
        else "Close abstract_compact_source_spectral_projection_theorem.py.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract compact source projection convergence terminal theorem")
    print(f"  projection convergence closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
