#!/usr/bin/env python3
r"""Terminal consequence for endpoint Riccati/Krawczyk capacities.

The endpoint coefficient enclosure only needs the finite Krawczyk capacity
numbers that bound coefficient and boundary-row perturbations.  This wrapper
exports those numbers without pulling the full Riccati/collocation proof
ledger into upstream audits.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def positive_number(data: dict[str, Any], key: str) -> bool:
    try:
        return float(data[key]) > 0.0
    except Exception:
        return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--krawczyk-json",
        default="endpoint_riccati_krawczyk_collocation.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_riccati_krawczyk_consequence_theorem.json",
    )
    args = parser.parse_args()

    kraw = load(args.krawczyk_json)
    coeff_ok = positive_number(kraw, "coefficientRadiusCapacityScaledCompanion")
    boundary_ok = positive_number(kraw, "boundaryRadiusCapacityUnscaledRows")
    fraction_ok = positive_number(kraw, "simultaneousCapacityFraction")
    endpoint_capacity = kraw.get(
        "rawEndpointEntryCapacityForChartThreshold",
        kraw.get("endpointEntryRadiusCapacity"),
    )
    try:
        endpoint_ok = float(endpoint_capacity) > 0.0
    except Exception:
        endpoint_ok = False
    closed = coeff_ok and boundary_ok and fraction_ok and endpoint_ok

    data: dict[str, Any] = {
        "theoremName": "endpoint Riccati Krawczyk capacity consequence theorem",
        "proofClass": "symbolic identity",
        "krawczykJson": args.krawczyk_json,
        "coefficientRadiusCapacityScaledCompanion": kraw.get(
            "coefficientRadiusCapacityScaledCompanion"
        ),
        "boundaryRadiusCapacityUnscaledRows": kraw.get(
            "boundaryRadiusCapacityUnscaledRows"
        ),
        "simultaneousCapacityFraction": kraw.get("simultaneousCapacityFraction"),
        "endpointEntryRadiusCapacity": endpoint_capacity,
        "actualKrawczykQ": kraw.get("actualKrawczykQ"),
        "endpointRiccatiKrawczykCapacityConsequenceClosed": closed,
        "statuses": {
            "coefficientCapacityStatus": status(
                "positive scaled companion coefficient capacity",
                coeff_ok,
                "The Riccati/Krawczyk collocation certificate reports a positive coefficient capacity.",
            ),
            "boundaryCapacityStatus": status(
                "positive active boundary-row capacity",
                boundary_ok,
                "The Riccati/Krawczyk collocation certificate reports a positive boundary-row capacity.",
            ),
            "simultaneousCapacityStatus": status(
                "positive simultaneous Krawczyk capacity fraction",
                fraction_ok,
                "The Riccati/Krawczyk collocation certificate reports a positive simultaneous perturbation fraction.",
            ),
            "endpointChartCapacityStatus": status(
                "positive endpoint-entry chart capacity",
                endpoint_ok,
                "The Riccati/Krawczyk collocation certificate reports a positive endpoint-entry chart radius.",
            ),
        },
        "proof": [
            "Import the finite endpoint Riccati/Krawczyk collocation certificate.",
            "Read off the positive coefficient, boundary, simultaneous, and endpoint-entry capacities.",
            "Export only those capacities to the endpoint coefficient interval enclosure.",
        ],
        "remainingAnalyticGap": None
        if closed
        else "One endpoint Krawczyk capacity is missing or nonpositive.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint Riccati Krawczyk consequence theorem")
    print(f"  coefficient capacity: {coeff_ok}")
    print(f"  boundary capacity: {boundary_ok}")
    print(f"  simultaneous capacity: {fraction_ok}")
    print(f"  endpoint chart capacity: {endpoint_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
