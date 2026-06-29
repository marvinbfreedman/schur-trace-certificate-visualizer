#!/usr/bin/env python3
r"""Normalized exterior/Grassmannian diagnostic for the endpoint rank shortcut.

The raw endpoint map

    M = R_R(s0) - R_L(s0)

is a 2 x 8 matrix.  Its full active row rank is equivalent to nonvanishing of
the Pluecker vector

    p_ij = det M[:, (i,j)],        0 <= i < j < 8.

The previous row-flow diagnostic tracked raw minors, whose scale is still
unstable.  This script moves to the projective Grassmannian coordinates:

    p_hat = p / ||p||_2.

This separates the question "does the active two-plane collapse?" from the
ill-conditioned absolute endpoint-flow scale.  It is still a numerical center,
not an interval proof, but it identifies the affine chart in which a future
ball certificate should be made.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def det2(a: float, b: float, c: float, d: float) -> float:
    return a * d - b * c


def plucker_from_matrix(mat: list[list[float]]) -> dict:
    if len(mat) != 2:
        raise ValueError("endpoint map must have two active rows")
    cols = len(mat[0])
    coordinates = []
    values = []
    for i, j in itertools.combinations(range(cols), 2):
        det = det2(mat[0][i], mat[0][j], mat[1][i], mat[1][j])
        values.append(det)
        coordinates.append(
            {
                "columns": [i, j],
                "determinant": det,
                "absDeterminant": abs(det),
            }
        )
    norm = math.sqrt(sum(value * value for value in values))
    if norm == 0.0:
        normalized = [0.0 for _ in values]
    else:
        normalized = [value / norm for value in values]
    for item, value in zip(coordinates, normalized):
        item["normalizedCoordinate"] = value
        item["absNormalizedCoordinate"] = abs(value)
    best = max(coordinates, key=lambda item: item["absNormalizedCoordinate"])
    return {
        "wedgeNorm": norm,
        "coordinates": coordinates,
        "normalizedVector": normalized,
        "bestNormalizedMinor": best,
        "topNormalizedMinors": sorted(
            coordinates,
            key=lambda item: item["absNormalizedCoordinate"],
            reverse=True,
        )[:8],
    }


def choose_persistent_chart(pluckers: list[dict]) -> dict:
    if not pluckers:
        raise ValueError("no Pluecker data")
    count = len(pluckers[0]["coordinates"])
    best_index = None
    best_margin = -1.0
    best_mean = -1.0
    for idx in range(count):
        values = [
            abs(plucker["coordinates"][idx]["normalizedCoordinate"])
            for plucker in pluckers
        ]
        margin = min(values)
        mean = sum(values) / len(values)
        if margin > best_margin or (margin == best_margin and mean > best_mean):
            best_index = idx
            best_margin = margin
            best_mean = mean
    item = pluckers[0]["coordinates"][best_index]
    return {
        "columns": item["columns"],
        "coordinateIndex": best_index,
        "minAbsNormalizedCoordinate": best_margin,
        "meanAbsNormalizedCoordinate": best_mean,
    }


def oriented_vectors(pluckers: list[dict]) -> list[list[float]]:
    out = []
    previous = None
    for plucker in pluckers:
        vector = list(plucker["normalizedVector"])
        if previous is not None:
            dot = sum(a * b for a, b in zip(previous, vector))
            if dot < 0:
                vector = [-value for value in vector]
        out.append(vector)
        previous = vector
    return out


def chart_coordinates(plucker: dict, chart_index: int) -> list[dict]:
    pivot = plucker["coordinates"][chart_index]["normalizedCoordinate"]
    out = []
    for idx, item in enumerate(plucker["coordinates"]):
        value = item["normalizedCoordinate"]
        out.append(
            {
                "columns": item["columns"],
                "coordinateIndex": idx,
                "ratioToChartMinor": value / pivot if pivot else None,
            }
        )
    return out


def compare_projective(left: dict, right: dict, left_vec: list[float], right_vec: list[float]) -> dict:
    dot = sum(a * b for a, b in zip(left_vec, right_vec))
    dot = max(-1.0, min(1.0, dot))
    abs_dot = abs(dot)
    projective_chordal = math.sqrt(max(0.0, 1.0 - abs_dot * abs_dot))
    oriented_chordal = math.sqrt(
        sum((a - b) * (a - b) for a, b in zip(left_vec, right_vec))
    )
    return {
        "leftOrder": left["chebyshevOrder"],
        "rightOrder": right["chebyshevOrder"],
        "orientedDot": dot,
        "projectiveAbsDot": abs_dot,
        "projectiveChordalDistance": projective_chordal,
        "orientedChordalDistance": oriented_chordal,
    }


def build_report(args: argparse.Namespace) -> dict:
    row_flow = load_json(args.row_flow_json)
    centers = row_flow.get("centers", [])
    pluckers = [plucker_from_matrix(center["endpointMap"]) for center in centers]
    chart = choose_persistent_chart(pluckers)
    oriented = oriented_vectors(pluckers)
    enriched = []
    for center, plucker, vector in zip(centers, pluckers, oriented):
        chart_values = chart_coordinates(plucker, chart["coordinateIndex"])
        chart_norm = math.sqrt(
            sum(
                (item["ratioToChartMinor"] or 0.0) ** 2
                for item in chart_values
            )
        )
        enriched.append(
            {
                "chebyshevOrder": center["chebyshevOrder"],
                "endpointMapRank": center["endpointMapRank"],
                "sigmaMin": center.get("sigmaMin"),
                "endpointMapFrobenius": center.get("endpointMapFrobenius"),
                "wedgeNorm": plucker["wedgeNorm"],
                "bestNormalizedMinor": plucker["bestNormalizedMinor"],
                "topNormalizedMinors": plucker["topNormalizedMinors"],
                "chartMinor": plucker["coordinates"][chart["coordinateIndex"]],
                "chartCoordinateNorm": chart_norm,
                "orientedNormalizedVector": vector,
                "chartCoordinates": chart_values,
            }
        )
    comparisons = [
        compare_projective(left, right, left_vec, right_vec)
        for left, right, left_vec, right_vec in zip(
            centers[:-1],
            centers[1:],
            oriented[:-1],
            oriented[1:],
        )
    ]
    tail_comparisons = comparisons[-args.tail_pairs :] if comparisons else []
    chart_margin = chart["minAbsNormalizedCoordinate"]
    tail_distance = max(
        (row["projectiveChordalDistance"] for row in tail_comparisons),
        default=float("inf"),
    )
    candidate_closed = bool(chart_margin >= args.chart_margin)
    tail_stable = bool(tail_comparisons and tail_distance <= args.projective_tol)
    return {
        "theoremName": "endpoint normalized exterior Grassmann center",
        "sourceRowFlow": "endpoint adjoint row-flow center",
        "orders": [center["chebyshevOrder"] for center in centers],
        "persistentChart": chart,
        "chartMarginThreshold": args.chart_margin,
        "projectiveTailTolerance": args.projective_tol,
        "tailPairsUsed": args.tail_pairs,
        "grassmannCenters": enriched,
        "projectiveComparisons": comparisons,
        "persistentChartNoncollapseStatus": status(
            "persistent affine Grassmann chart identified",
            candidate_closed,
            (
                "Closed at the numerical-center level if one Pluecker "
                "coordinate stays uniformly away from zero across all recorded "
                "orders."
            ),
        ),
        "tailProjectiveStabilityStatus": status(
            "tail projective Grassmann stability",
            tail_stable,
            (
                "Closed at the numerical-center level if the latest normalized "
                "two-plane comparisons are below the configured projective "
                "distance tolerance."
            ),
        ),
        "rankShortcutReadyStatus": status(
            "full-rank shortcut ready for interval proof",
            False,
            (
                "Still open.  This diagnostic identifies the projective chart "
                "for a future ball certificate, but it does not enclose the "
                "exact endpoint flow."
            ),
        ),
        "conclusion": (
            "The raw endpoint-map scale is unstable, but the normalized "
            "Pluecker coordinate in the persistent chart is large.  The "
            "latest projective comparison is the relevant convergence metric "
            "for the full-rank shortcut."
        ),
        "nextTarget": (
            "Build an interval enclosure directly in the persistent Pluecker "
            "chart, or integrate the corresponding Riccati chart equations "
            "with ball arithmetic."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--row-flow-json", default="endpoint_adjoint_row_flow_center.json")
    parser.add_argument("--json-out", default="endpoint_grassmann_flow_center.json")
    parser.add_argument("--chart-margin", type=float, default=0.5)
    parser.add_argument("--projective-tol", type=float, default=0.05)
    parser.add_argument("--tail-pairs", type=int, default=1)
    args = parser.parse_args()

    data = build_report(args)
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint normalized exterior Grassmann center")
    print(f"  persistent chart: {data['persistentChart']['columns']}")
    print(
        "  chart min |p_hat|: "
        f"{data['persistentChart']['minAbsNormalizedCoordinate']:.12e}"
    )
    for center in data["grassmannCenters"]:
        best = center["bestNormalizedMinor"]
        print(
            f"  order={center['chebyshevOrder']} "
            f"best={best['columns']} "
            f"|p_hat|={best['absNormalizedCoordinate']:.12e} "
            f"wedge_norm={center['wedgeNorm']:.12e}"
        )
    for row in data["projectiveComparisons"]:
        print(
            f"  {row['leftOrder']}->{row['rightOrder']} "
            f"abs_dot={row['projectiveAbsDot']:.12e} "
            f"dist={row['projectiveChordalDistance']:.12e}"
        )
    print(
        "  persistent chart noncollapse: "
        f"{data['persistentChartNoncollapseStatus']['closed']}"
    )
    print(
        "  tail projective stable: "
        f"{data['tailProjectiveStabilityStatus']['closed']}"
    )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
