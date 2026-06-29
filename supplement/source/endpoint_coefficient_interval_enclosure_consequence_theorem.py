#!/usr/bin/env python3
r"""Terminal consequence of the endpoint coefficient interval enclosure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {"label": label, "closed": closed, "status": "closed" if closed else "open", "reason": reason}
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coefficient-json",
        default="endpoint_coefficient_interval_enclosure.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_coefficient_interval_enclosure_consequence_theorem.json",
    )
    args = parser.parse_args()

    coeff = load(args.coefficient_json)
    scaled = bool(coeff.get("scaledCompanionIntervalStatus", {}).get("closed"))
    boundary = bool(coeff.get("boundaryRowIntervalStatus", {}).get("closed"))
    budget = bool(coeff.get("simultaneousKrawczykInputStatus", {}).get("closed"))
    shortcut = bool(coeff.get("fullRankShortcutStatus", {}).get("closed"))
    closed = scaled and boundary and budget and shortcut

    data = {
        "theoremName": "endpoint coefficient interval enclosure consequence theorem",
        "proofClass": "symbolic identity",
        "scaledCompanionIntervalStatus": status(
            "scaled companion interval enclosure",
            scaled,
            "The scaled companion coefficients are enclosed below the finite Krawczyk capacity.",
        ),
        "boundaryRowIntervalStatus": status(
            "boundary row interval enclosure",
            boundary,
            "The active endpoint boundary rows are enclosed below the finite Krawczyk capacity.",
        ),
        "simultaneousKrawczykInputStatus": status(
            "simultaneous Krawczyk input",
            budget,
            "The simultaneous coefficient and boundary-row interval budgets close.",
        ),
        "fullRankShortcutStatus": status(
            "full-rank shortcut",
            shortcut,
            "The interval coefficient enclosures fit inside the full-rank shortcut capacity.",
        ),
        "scaledCompanionIntervalRadius": coeff.get("scaledCompanionIntervalRadius"),
        "boundaryRowIntervalRadius": coeff.get("boundaryRowIntervalRadius"),
        "scaledCompanionCapacity": coeff.get("scaledCompanionCapacity"),
        "boundaryRowCapacity": coeff.get("boundaryRowCapacity"),
        "coefficientRadiusScale": coeff.get("scaledCompanionIntervalRadius"),
        "endpointCoefficientIntervalEnclosureClosed": closed,
        "remainingAnalyticGap": None if closed else "Close endpoint_coefficient_interval_enclosure.py.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint coefficient interval enclosure consequence theorem")
    print(f"  coefficient enclosure closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
