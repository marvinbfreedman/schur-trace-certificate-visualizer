#!/usr/bin/env python3
r"""Ball certificate in the persistent Pluecker chart.

The normalized exterior diagnostic identifies the affine chart p_01 != 0.
This script turns that into an explicit ball statement:

    if min(||x-v||_2, ||x+v||_2) <= r

where v is the normalized Pluecker center and x is the exact normalized
Pluecker vector, then

    |x_01| >= |v_01| - r.

Thus any verified projective ball with r < |v_01| proves full rank, and any
ball with r <= |v_01|-1/2 proves the stronger margin |x_01| >= 1/2.

This file does not yet enclose the exact endpoint ODE/Riccati flow.  It is the
chart-ball algebra and the quantitative target for that enclosure.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def sci(x: float) -> float:
    return float(x)


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def iv_abs_sup(x) -> mp.mpf:
    return max(abs(mp.mpf(x.a)), abs(mp.mpf(x.b)))


def iv_abs_inf(x) -> mp.mpf:
    lo = mp.mpf(x.a)
    hi = mp.mpf(x.b)
    if lo <= 0 <= hi:
        return mp.mpf("0")
    return min(abs(lo), abs(hi))


def iv_ball(x: float, radius: mp.mpf):
    x = mp.mpf(str(x))
    return mp.iv.mpf([x - radius, x + radius])


def raw_entry_ratio_lower(mat: list[list[float]], chart_cols: list[int], radius: mp.mpf) -> mp.mpf:
    minors = []
    pivot = None
    cols = len(mat[0])
    for i in range(cols):
        for j in range(i + 1, cols):
            a = iv_ball(mat[0][i], radius)
            b = iv_ball(mat[0][j], radius)
            c = iv_ball(mat[1][i], radius)
            d = iv_ball(mat[1][j], radius)
            det = a * d - b * c
            minors.append(det)
            if [i, j] == chart_cols:
                pivot = det
    if pivot is None:
        raise ValueError(f"chart columns not found: {chart_cols}")
    numerator = iv_abs_inf(pivot)
    denominator = mp.sqrt(mp.fsum(iv_abs_sup(det) ** 2 for det in minors))
    if denominator == 0:
        return mp.mpf("0")
    return numerator / denominator


def raw_entry_radius_capacity(
    mat: list[list[float]],
    chart_cols: list[int],
    threshold: mp.mpf,
) -> dict:
    hi = mp.mpf("1")
    while raw_entry_ratio_lower(mat, chart_cols, hi) >= threshold:
        hi *= 2
        if hi > mp.mpf("1e12"):
            break
    lo = mp.mpf("0")
    for _ in range(80):
        mid = (lo + hi) / 2
        if raw_entry_ratio_lower(mat, chart_cols, mid) >= threshold:
            lo = mid
        else:
            hi = mid
    max_entry = max(abs(value) for row in mat for value in row)
    return {
        "threshold": sci(threshold),
        "uniformAbsoluteEntryRadius": sci(lo),
        "relativeToMaxEntry": sci(lo / mp.mpf(str(max_entry))) if max_entry else None,
        "ratioLowerAtRadius": sci(raw_entry_ratio_lower(mat, chart_cols, lo)),
    }


def orient_to_positive_chart(vector: list[float], chart_index: int) -> list[float]:
    if vector[chart_index] < 0:
        return [-value for value in vector]
    return list(vector)


def chart_ball(center: list[float], chart_index: int, radius: float, threshold: float) -> dict:
    pivot = abs(center[chart_index])
    lower = pivot - radius
    upper = min(1.0, pivot + radius)
    noncollapse = lower > 0
    threshold_closed = lower >= threshold
    ratios = []
    for idx, value in enumerate(center):
        center_ratio = value / center[chart_index]
        if idx == chart_index:
            diff_bound = 0.0
            abs_upper = 1.0
        elif lower > 0:
            diff_bound = radius * (pivot + abs(value)) / (pivot * lower)
            abs_upper = (abs(value) + radius) / lower
        else:
            diff_bound = math.inf
            abs_upper = math.inf
        ratios.append(
            {
                "coordinateIndex": idx,
                "centerRatio": center_ratio,
                "ratioDiffBound": diff_bound,
                "ratioInterval": [
                    center_ratio - diff_bound,
                    center_ratio + diff_bound,
                ] if math.isfinite(diff_bound) else None,
                "absRatioUpper": abs_upper,
            }
        )
    return {
        "projectiveRadius": radius,
        "pivotAbsCenter": pivot,
        "pivotAbsLower": lower,
        "pivotAbsUpper": upper,
        "provesNoncollapse": noncollapse,
        "provesThreshold": threshold_closed,
        "threshold": threshold,
        "maxChartRatioDiffBound": max(
            row["ratioDiffBound"]
            for row in ratios
            if math.isfinite(row["ratioDiffBound"])
        ),
        "maxAbsChartRatioUpper": max(
            row["absRatioUpper"]
            for row in ratios
            if math.isfinite(row["absRatioUpper"])
        ),
        "chartRatioBounds": ratios,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grassmann-json", default="endpoint_grassmann_center_consequence_theorem.json")
    parser.add_argument("--row-flow-json", default="")
    parser.add_argument("--center-order", type=int, default=0)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--safety-factors", default="1 2 4 8 12 16")
    parser.add_argument("--cert-radius", type=float, default=0.125)
    parser.add_argument("--iv-dps", type=int, default=70)
    parser.add_argument("--json-out", default="endpoint_grassmann_chart_ball_certificate.json")
    args = parser.parse_args()

    mp.iv.dps = args.iv_dps
    data = load_json(args.grassmann_json)
    centers = data["grassmannCenters"]
    center = centers[-1]
    if args.center_order:
        matches = [row for row in centers if row["chebyshevOrder"] == args.center_order]
        if not matches:
            raise SystemExit(f"center order not found: {args.center_order}")
        center = matches[0]

    chart = data["persistentChart"]
    chart_index = int(chart["coordinateIndex"])
    chart_cols = list(chart["columns"])
    vector = orient_to_positive_chart(center["orientedNormalizedVector"], chart_index)
    pivot = abs(vector[chart_index])
    threshold = float(args.threshold)
    radius_capacity_noncollapse = pivot
    radius_capacity_threshold = max(0.0, pivot - threshold)
    comparisons = data.get("projectiveComparisons", [])
    latest_distance = float(comparisons[-1]["projectiveChordalDistance"]) if comparisons else math.inf
    max_tail_distance = max(
        (float(row["projectiveChordalDistance"]) for row in comparisons),
        default=math.inf,
    )
    factors = [float(piece) for piece in args.safety_factors.replace(",", " ").split()]
    empirical_balls = [
        chart_ball(vector, chart_index, latest_distance * factor, threshold)
        for factor in factors
    ]
    certified_ball = chart_ball(vector, chart_index, args.cert_radius, threshold)

    row_flow_path = args.row_flow_json or data.get("sourceRowFlowJson")
    raw_capacity = data.get("precomputedRawEntryBallCapacity")
    if row_flow_path:
        row_flow = load_json(row_flow_path)
        raw_centers = row_flow.get("centers", [])
        raw_center = raw_centers[-1]
        if args.center_order:
            matches = [
                row for row in raw_centers
                if row["chebyshevOrder"] == args.center_order
            ]
            if matches:
                raw_center = matches[0]
        raw_capacity = raw_entry_radius_capacity(
            raw_center["endpointMap"],
            chart_cols,
            mp.mpf(str(threshold)),
        )

    exact_radius_ready = certified_ball["provesThreshold"]
    report = {
        "theoremName": "endpoint Grassmann chart ball certificate",
        "sourceGrassmannJson": args.grassmann_json,
        "sourceRowFlowJson": row_flow_path,
        "centerOrder": center["chebyshevOrder"],
        "persistentChart": chart,
        "chartCoordinate": chart_cols,
        "pivotAbsCenter": pivot,
        "threshold": threshold,
        "radiusCapacityForNoncollapse": radius_capacity_noncollapse,
        "radiusCapacityForThreshold": radius_capacity_threshold,
        "latestProjectiveDistance": latest_distance,
        "maxRecordedProjectiveDistance": max_tail_distance,
        "certifiedProjectiveRadius": args.cert_radius,
        "certifiedProjectiveBall": certified_ball,
        "empiricalSafetyBalls": empirical_balls,
        "rawEntryBallCapacity": raw_capacity,
        "conditionalChartBallTheoremStatus": status(
            "conditional chart-ball noncollapse theorem",
            exact_radius_ready,
            (
                "If an exact interval/Riccati computation encloses the "
                "normalized Pluecker vector inside the certified projective "
                "ball, then p_01 is nonzero and the stronger threshold margin "
                "holds."
            ),
        ),
        "exactFlowEnclosureStatus": status(
            "exact endpoint Riccati/Pluecker flow enclosure",
            False,
            (
                "Open.  This script supplies the chart-ball target and interval "
                "algebra, but it does not yet enclose the exact ODE flow."
            ),
        ),
        "conclusion": (
            "The persistent chart has large margin.  A projective enclosure "
            f"radius below {radius_capacity_threshold:.6g} proves "
            f"|p_hat_01| >= {threshold:.6g}; the latest observed projective "
            f"distance is {latest_distance:.6g}."
        ),
        "nextTarget": (
            "Integrate the Riccati equations for the chart ratios p_J/p_01 "
            "with ball arithmetic and prove the exact projective radius is "
            "below the certified radius."
        ),
    }
    Path(args.json_out).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("Endpoint Grassmann chart ball certificate")
    print(f"  center order: {center['chebyshevOrder']}")
    print(f"  chart: {chart_cols}")
    print(f"  |p_hat_chart| center: {pivot:.12e}")
    print(f"  radius for noncollapse: {radius_capacity_noncollapse:.12e}")
    print(f"  radius for |p_hat_chart| >= {threshold}: {radius_capacity_threshold:.12e}")
    print(f"  latest projective distance: {latest_distance:.12e}")
    print(f"  certified radius: {args.cert_radius:.12e}")
    print(f"  certified lower: {certified_ball['pivotAbsLower']:.12e}")
    if raw_capacity:
        print(
            "  raw entry radius for same threshold: "
            f"{raw_capacity['uniformAbsoluteEntryRadius']:.12e}"
        )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
