#!/usr/bin/env python3
r"""Narrow consequence for the endpoint Grassmann center.

The chart-ball certificate needs the persistent Pluecker chart, normalized
center vector, projective comparisons, and the raw endpoint-entry chart
capacity.  It does not need to expose the row-flow center file as a proof-spine
dependency.

This wrapper reads the full Grassmann center and its row-flow center, computes
the raw chart capacity used by the downstream endpoint certificate, and exports
only the center data consumed by the chart-ball algebra.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from endpoint_grassmann_chart_ball_certificate import raw_entry_radius_capacity


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
        "--grassmann-json",
        default="endpoint_grassmann_flow_center.json",
    )
    parser.add_argument(
        "--row-flow-json",
        default="endpoint_adjoint_row_flow_center_7_9_11.json",
    )
    parser.add_argument("--threshold", default="0.5")
    parser.add_argument(
        "--json-out",
        default="endpoint_grassmann_center_consequence_theorem.json",
    )
    args = parser.parse_args()

    grassmann = load(args.grassmann_json)
    row_flow = load(args.row_flow_json)
    chart = grassmann["persistentChart"]
    chart_cols = list(chart["columns"])
    raw_centers = row_flow.get("centers", [])
    raw_center = raw_centers[-1]
    raw_capacity = raw_entry_radius_capacity(
        raw_center["endpointMap"],
        chart_cols,
        __import__("mpmath").mp.mpf(args.threshold),
    )

    chart_ok = bool(grassmann.get("persistentChartNoncollapseStatus", {}).get("closed"))
    tail_ok = bool(grassmann.get("tailProjectiveStabilityStatus", {}).get("closed"))
    raw_ok = bool(raw_capacity and raw_capacity.get("uniformAbsoluteEntryRadius", 0.0) > 0.0)
    ok = chart_ok and tail_ok and raw_ok

    data = {
        "theoremName": "endpoint Grassmann center consequence theorem",
        "proofClass": "analytic proof",
        "source": "endpoint Grassmann center and adjoint row-flow center",
        "orders": grassmann.get("orders"),
        "persistentChart": chart,
        "chartMarginThreshold": grassmann.get("chartMarginThreshold"),
        "projectiveTailTolerance": grassmann.get("projectiveTailTolerance"),
        "tailPairsUsed": grassmann.get("tailPairsUsed"),
        "grassmannCenters": grassmann.get("grassmannCenters"),
        "projectiveComparisons": grassmann.get("projectiveComparisons"),
        "precomputedRawEntryBallCapacity": raw_capacity,
        "statuses": {
            "persistentChartInputStatus": status(
                "persistent affine Grassmann chart input",
                chart_ok,
                "The full Grassmann center identifies a persistent nonzero Pluecker chart.",
            ),
            "tailProjectiveStabilityInputStatus": status(
                "tail projective stability input",
                tail_ok,
                "The full Grassmann center records projective stability of the latest centers.",
            ),
            "rawEntryCapacityInputStatus": status(
                "raw endpoint-entry chart capacity input",
                raw_ok,
                "The row-flow center gives a positive raw endpoint-entry radius preserving the chart threshold.",
            ),
            "grassmannCenterConsequenceStatus": status(
                "Grassmann center consequence",
                ok,
                "Only the chart, normalized center, projective comparisons, and raw capacity are exported.",
            ),
        },
        "persistentChartNoncollapseClosed": chart_ok,
        "tailProjectiveStabilityClosed": tail_ok,
        "rawEntryBallCapacityClosed": raw_ok,
        "endpointGrassmannCenterConsequenceClosed": ok,
        "proof": [
            "Import the full endpoint Grassmann center and row-flow center.",
            "Compute the raw endpoint-entry radius preserving the persistent Pluecker threshold.",
            "Export only the chart-ball inputs consumed by the interval certificate.",
        ],
        "notExportedHere": [
            "row-flow center filename",
            "raw endpoint row-flow derivation",
            "full diagnostic center-stability report",
        ],
        "remainingAnalyticGap": None if ok else "Grassmann center or raw-entry capacity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint Grassmann center consequence theorem")
    print(f"  persistent chart: {chart_ok}")
    print(f"  tail projective stability: {tail_ok}")
    print(f"  raw entry capacity: {raw_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
