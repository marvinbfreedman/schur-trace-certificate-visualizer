#!/usr/bin/env python3
r"""Terminal consequence of the endpoint Grassmann chart ball certificate."""

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
        "--chart-ball-json",
        default="endpoint_grassmann_chart_ball_certificate.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_grassmann_chart_ball_consequence_theorem.json",
    )
    args = parser.parse_args()

    chart = load(args.chart_ball_json)
    conditional = bool(chart.get("conditionalChartBallTheoremStatus", {}).get("closed"))
    has_capacity = (
        chart.get("chartCoordinate") is not None
        and chart.get("threshold") is not None
        and chart.get("radiusCapacityForNoncollapse") is not None
    )
    closed = conditional and has_capacity

    data = {
        "theoremName": "endpoint Grassmann chart ball consequence theorem",
        "proofClass": "symbolic identity",
        "chartCapacityInputStatus": status(
            "Grassmann chart capacity input",
            has_capacity,
            (
                "The chart coordinate, threshold, and noncollapse capacity "
                "needed by the synchronized endpoint coefficient certificate "
                "are available."
            ),
        ),
        "conditionalChartBallTheoremStatus": status(
            "conditional chart ball theorem",
            conditional,
            "The certified chart ball remains inside the noncollapse capacity.",
        ),
        "chartCoordinate": chart.get("chartCoordinate"),
        "threshold": chart.get("threshold"),
        "radiusCapacityForNoncollapse": chart.get("radiusCapacityForNoncollapse"),
        "radiusCapacityForThreshold": chart.get("radiusCapacityForThreshold"),
        "certifiedProjectiveRadius": chart.get("certifiedProjectiveRadius"),
        "endpointGrassmannChartCapacityInputClosed": closed,
        "remainingAnalyticGap": None
        if closed
        else "Close the conditional chart-ball capacity input.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint Grassmann chart ball consequence theorem")
    print(f"  chart ball closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
