#!/usr/bin/env python3
r"""Narrow active-projection convergence consequence theorem.

This wrapper exposes only the compact-source spectral consequence needed by
the abstract min-max tail passage:

    active Riesz projections converge in norm,
    hence inactive projections also converge in norm.

The compact compression and Riesz-contour proof remains inside
``abstract_compact_source_spectral_projection_theorem.json``.
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
        default="abstract_compact_source_projection_convergence_terminal_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_active_projection_convergence_consequence_theorem.json",
    )
    args = parser.parse_args()

    spectral = load(args.spectral_json)
    riesz_projection_ok = bool(
        spectral.get("compactSourceRieszProjectionConvergenceClosed")
        or nested_closed(spectral, "statuses", "rieszProjectionStabilityStatus")
    )
    splitting_ok = bool(
        spectral.get("sourceActiveInactiveSplittingStabilityClosed")
        or nested_closed(spectral, "statuses", "activeInactiveSplittingStabilityStatus")
    )
    active_projection_ok = riesz_projection_ok and splitting_ok
    inactive_projection_ok = active_projection_ok

    data = {
        "theoremName": "abstract active projection convergence consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The compact-source spectral projection theorem implies norm "
            "convergence of active Riesz projections and, by subtraction from "
            "the identity, norm convergence of inactive projections."
        ),
        "statuses": {
            "activeRieszProjectionConvergenceStatus": status(
                "active Riesz projection convergence",
                active_projection_ok,
                (
                    "Imported from the compact-source spectral projection "
                    "theorem: the isolated active source cluster has stable "
                    "Riesz projections and active/inactive splitting."
                ),
            ),
            "inactiveProjectionConvergenceStatus": status(
                "inactive projection convergence",
                inactive_projection_ok,
                (
                    "If P_active,N -> P_active in norm, then "
                    "I-P_active,N -> I-P_active in norm with the same norm "
                    "error."
                ),
            ),
        },
        "activeRieszProjectionConvergenceClosed": active_projection_ok,
        "inactiveProjectionConvergenceClosed": inactive_projection_ok,
        "compactSourceRieszProjectionConvergenceClosed": active_projection_ok,
        "sourceActiveInactiveSplittingStabilityClosed": splitting_ok,
        "proof": [
            "Import compact-source Riesz projection convergence.",
            "Import stability of the active/inactive source splitting.",
            "Subtract active projections from the identity to get inactive projection convergence.",
            "Expose only the projection-convergence consequence needed by the min-max tail passage.",
        ],
        "remainingAnalyticGap": None
        if inactive_projection_ok
        else "Compact-source active projection convergence input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract active projection convergence consequence theorem")
    print(f"  active projection convergence: {active_projection_ok}")
    print(f"  inactive projection convergence: {inactive_projection_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
