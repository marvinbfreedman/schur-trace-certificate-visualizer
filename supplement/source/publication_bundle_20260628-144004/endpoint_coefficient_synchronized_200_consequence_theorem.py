#!/usr/bin/env python3
r"""Terminal consequence of the synchronized 200-point endpoint certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def item_closed(data: dict, key: str) -> bool:
    return bool(data.get(key, {}).get("closed"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {"label": label, "closed": closed, "status": "closed" if closed else "open", "reason": reason}
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--endpoint-json",
        default="endpoint_coefficient_synchronized_200_certificate.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_coefficient_synchronized_200_consequence_theorem.json",
    )
    args = parser.parse_args()

    endpoint = load(args.endpoint_json)
    center = item_closed(endpoint, "centerSynchronizationStatus")
    companion = item_closed(endpoint, "analyticCompanionInputStatus")
    boundary = item_closed(endpoint, "analyticBoundaryInputStatus")
    krawczyk = item_closed(endpoint, "synchronizedKrawczykInputStatus")
    q = float(endpoint.get("actualKrawczykQ", 1.0))
    chart = float(endpoint.get("endpointRadiusOverChartCapacity", 1.0))
    closed = center and companion and boundary and krawczyk and q < 1.0 and chart < 1.0

    data = {
        "theoremName": "synchronized 200-point endpoint coefficient consequence theorem",
        "proofClass": "symbolic identity",
        "centerSynchronizationStatus": status(
            "synchronized 200-point center",
            center,
            "The synchronized 200-point endpoint center is the certified center used by downstream interval inputs.",
        ),
        "analyticCompanionInputStatus": status(
            "analytic companion coefficient interval input",
            companion,
            "The scaled companion coefficients are enclosed inside the Krawczyk input budget.",
        ),
        "analyticBoundaryInputStatus": status(
            "analytic active boundary row interval input",
            boundary,
            "The active endpoint boundary rows are enclosed inside the Krawczyk input budget.",
        ),
        "synchronizedKrawczykInputStatus": status(
            "synchronized Krawczyk input",
            krawczyk,
            "The synchronized interval inputs close the persistent Grassmann/Krawczyk chart certificate.",
        ),
        "actualKrawczykQ": endpoint.get("actualKrawczykQ"),
        "endpointRadiusOverChartCapacity": endpoint.get("endpointRadiusOverChartCapacity"),
        "scaledCompanionRadiusOverSimultaneousBudget": endpoint.get(
            "scaledCompanionRadiusOverSimultaneousBudget"
        ),
        "boundaryRowRadiusOverSimultaneousBudget": endpoint.get(
            "boundaryRowRadiusOverSimultaneousBudget"
        ),
        "endpointCoefficientSynchronized200Closed": closed,
        "remainingAnalyticGap": None if closed else "Close endpoint_coefficient_synchronized_200_certificate.py.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint coefficient synchronized 200 consequence theorem")
    print(f"  endpoint coefficient input closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
