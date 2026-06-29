#!/usr/bin/env python3
r"""Stability ledger for endpoint-map centers.

The rank-ball theorem can prove endpoint compatibility only after we choose a
trustworthy numerical center M0 for the exact endpoint map.  This ledger
compares several finite constructions:

* sampled moving eigenrow + finite-difference coefficient field,
* exact local confluent eigen-derivatives on each ODE segment.

Large variation between centers means the desired inequality

    ||M_exact-M0||_F < sigma_min(M0)

has not yet been reduced to a rigorous enclosure problem; the center itself
must first be rebuilt using a stable/validated ODE integration scheme.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    candidate = Path(path)
    if not candidate.exists():
        return None
    return json.loads(candidate.read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def summarize(path: str) -> dict:
    data = load_optional(path)
    if not data:
        return {"path": path, "present": False}
    eigs = data.get("endpointMapGramEigenvalues") or []
    sigmas = [max(0.0, float(value)) ** 0.5 for value in eigs]
    sigma_min = min(sigmas) if sigmas else None
    sigma_max = max(sigmas) if sigmas else None
    return {
        "path": path,
        "present": True,
        "theoremName": data.get("theoremName"),
        "odeSamples": data.get("odeSamples"),
        "activeDim": data.get("activeDim"),
        "endpointMapRank": data.get("endpointMapRank"),
        "endpointMapFrobenius": data.get("endpointMapFrobenius"),
        "flowLeftFrobenius": data.get("flowLeftFrobenius"),
        "flowRightFrobenius": data.get("flowRightFrobenius"),
        "maxEndpointResidualNorm": data.get("maxEndpointResidualNorm"),
        "maxSelectedZNorm": data.get("maxSelectedZNorm"),
        "sigmaMin": sigma_min,
        "sigmaMax": sigma_max,
        "conditionNumber": sigma_max / sigma_min if sigma_min and sigma_min > 0 else None,
    }


def ratio(a, b):
    if a is None or b is None or b == 0:
        return None
    return float(a) / float(b)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sampled-jsons",
        default=(
            "adjoint_green_endpoint_selection.json "
            "adjoint_green_endpoint_selection_o13.json "
            "adjoint_green_endpoint_selection_o15.json"
        ),
    )
    parser.add_argument(
        "--exact-segment-jsons",
        default=(
            "adjoint_green_endpoint_selection_exact_o5.json "
            "adjoint_green_endpoint_selection_exact_o7.json"
        ),
    )
    parser.add_argument(
        "--rank-json",
        default="endpoint_map_rank_ball_certificate_with_variants.json",
    )
    parser.add_argument(
        "--exact-rank-json",
        default="endpoint_map_rank_ball_certificate_exact_segments.json",
    )
    parser.add_argument(
        "--chebyshev-jsons",
        default="endpoint_flow_chebyshev_center.json endpoint_flow_chebyshev_center_7_9.json",
    )
    parser.add_argument(
        "--magnus-jsons",
        default="endpoint_flow_magnus_center_2_4.json",
    )
    parser.add_argument(
        "--row-flow-jsons",
        default=(
            "endpoint_adjoint_row_flow_center.json "
            "endpoint_adjoint_row_flow_center_7_9_11.json"
        ),
    )
    parser.add_argument(
        "--grassmann-jsons",
        default="endpoint_grassmann_flow_center.json",
    )
    parser.add_argument(
        "--chart-ball-json",
        default="endpoint_grassmann_chart_ball_certificate.json",
    )
    parser.add_argument(
        "--riccati-json",
        default="endpoint_riccati_flow_enclosure.json",
    )
    parser.add_argument(
        "--krawczyk-json",
        default="endpoint_riccati_krawczyk_collocation.json",
    )
    parser.add_argument(
        "--coefficient-ball-json",
        default="endpoint_coefficient_ball_enclosure.json",
    )
    parser.add_argument("--json-out", default="endpoint_map_center_stability.json")
    args = parser.parse_args()

    sampled_paths = args.sampled_jsons.replace(",", " ").split()
    exact_paths = args.exact_segment_jsons.replace(",", " ").split()
    sampled = [summarize(path) for path in sampled_paths]
    exact = [summarize(path) for path in exact_paths]
    rank = load_optional(args.rank_json) or {}
    exact_rank = load_optional(args.exact_rank_json) or {}
    cheb_paths = args.chebyshev_jsons.replace(",", " ").split()
    cheb_runs = [load_optional(path) for path in cheb_paths]
    magnus_paths = args.magnus_jsons.replace(",", " ").split()
    magnus_runs = [load_optional(path) for path in magnus_paths]
    row_flow_paths = args.row_flow_jsons.replace(",", " ").split()
    row_flow_runs = [load_optional(path) for path in row_flow_paths]
    grassmann_paths = args.grassmann_jsons.replace(",", " ").split()
    grassmann_runs = [load_optional(path) for path in grassmann_paths]
    chart_ball = load_optional(args.chart_ball_json) or {}
    riccati = load_optional(args.riccati_json) or {}
    krawczyk = load_optional(args.krawczyk_json) or {}
    coefficient_ball = load_optional(args.coefficient_ball_json) or {}
    cheb_summaries = []
    cheb_comparisons = []
    for path, run in zip(cheb_paths, cheb_runs):
        if not run:
            cheb_summaries.append({"path": path, "present": False})
            continue
        cheb_summaries.append(
            {
                "path": path,
                "present": True,
                "orders": run.get("orders"),
                "rankMargin": run.get("rankMargin"),
                "stabilityStatus": run.get("chebyshevCenterStabilityStatus", {}),
            }
        )
        for comparison in run.get("comparisons", []):
            item = dict(comparison)
            item["source"] = path
            cheb_comparisons.append(item)
    magnus_summaries = []
    magnus_comparisons = []
    for path, run in zip(magnus_paths, magnus_runs):
        if not run:
            magnus_summaries.append({"path": path, "present": False})
            continue
        magnus_summaries.append(
            {
                "path": path,
                "present": True,
                "steps": run.get("steps"),
                "rankMargin": run.get("rankMargin"),
                "stabilityStatus": run.get("magnusCenterStabilityStatus", {}),
            }
        )
        for comparison in run.get("comparisons", []):
            item = dict(comparison)
            item["source"] = path
            magnus_comparisons.append(item)
    row_flow_summaries = []
    row_flow_comparisons = []
    for path, run in zip(row_flow_paths, row_flow_runs):
        if not run:
            row_flow_summaries.append({"path": path, "present": False})
            continue
        row_flow_summaries.append(
            {
                "path": path,
                "present": True,
                "orders": run.get("orders"),
                "rankMargin": run.get("rankMargin"),
                "frobeniusStabilityStatus": run.get("frobeniusStabilityStatus", {}),
                "minorStabilityStatus": run.get("minorStabilityStatus", {}),
            }
        )
        for comparison in run.get("comparisons", []):
            item = dict(comparison)
            item["source"] = path
            row_flow_comparisons.append(item)
    grassmann_summaries = []
    grassmann_comparisons = []
    for path, run in zip(grassmann_paths, grassmann_runs):
        if not run:
            grassmann_summaries.append({"path": path, "present": False})
            continue
        grassmann_summaries.append(
            {
                "path": path,
                "present": True,
                "orders": run.get("orders"),
                "persistentChart": run.get("persistentChart"),
                "persistentChartNoncollapseStatus": run.get(
                    "persistentChartNoncollapseStatus",
                    {},
                ),
                "tailProjectiveStabilityStatus": run.get(
                    "tailProjectiveStabilityStatus",
                    {},
                ),
            }
        )
        for comparison in run.get("projectiveComparisons", []):
            item = dict(comparison)
            item["source"] = path
            grassmann_comparisons.append(item)

    sampled_frobs = [
        float(row["endpointMapFrobenius"])
        for row in sampled
        if row.get("present") and row.get("endpointMapFrobenius") is not None
    ]
    exact_frobs = [
        float(row["endpointMapFrobenius"])
        for row in exact
        if row.get("present") and row.get("endpointMapFrobenius") is not None
    ]
    sampled_spread = max(sampled_frobs) / min(sampled_frobs) if len(sampled_frobs) >= 2 else None
    exact_spread = max(exact_frobs) / min(exact_frobs) if len(exact_frobs) >= 2 else None

    sampled_stable = bool(sampled_spread is not None and sampled_spread < 2.0)
    exact_stable = bool(exact_spread is not None and exact_spread < 2.0)
    center_ready = sampled_stable or exact_stable
    cheb_stable = bool(
        cheb_comparisons
        and all(row.get("passesRankMargin") for row in cheb_comparisons)
    )
    magnus_stable = bool(
        magnus_comparisons
        and all(row.get("passesRankMargin") for row in magnus_comparisons)
    )
    row_flow_stable = bool(
        row_flow_comparisons
        and all(row.get("sameBestMinor") for row in row_flow_comparisons)
        and all(
            row.get("bestMinorRelativeDiff") is not None
            and float(row["bestMinorRelativeDiff"]) < 0.25
            for row in row_flow_comparisons
        )
    )
    grassmann_stable = bool(
        grassmann_summaries
        and all(
            row.get("persistentChartNoncollapseStatus", {}).get("closed")
            and row.get("tailProjectiveStabilityStatus", {}).get("closed")
            for row in grassmann_summaries
            if row.get("present")
        )
    )
    chart_ball_closed = bool(
        chart_ball
        and chart_ball.get("conditionalChartBallTheoremStatus", {}).get("closed")
    )
    riccati_equations_closed = bool(
        riccati
        and riccati.get("exteriorEquationStatus", {}).get("closed")
        and riccati.get("riccatiChartEquationStatus", {}).get("closed")
        and riccati.get("chartTubeNoncollapseStatus", {}).get("closed")
    )
    krawczyk_closed = bool(
        krawczyk
        and krawczyk.get("linearCollocationKrawczykStatus", {}).get("closed")
    )
    coefficient_ball_closed = bool(
        coefficient_ball
        and coefficient_ball.get("validatedKrawczykCoefficientInputStatus", {}).get("closed")
    )

    data = {
        "theoremName": "endpoint map center stability",
        "sampledFiniteDifferenceCenters": sampled,
        "exactSegmentDerivativeCenters": exact,
        "sampledFrobeniusSpread": sampled_spread,
        "exactSegmentFrobeniusSpread": exact_spread,
        "sampledVariantRankBall": rank,
        "exactSegmentVariantRankBall": exact_rank,
        "chebyshevCenters": cheb_summaries,
        "chebyshevComparisons": cheb_comparisons,
        "magnusCenters": magnus_summaries,
        "magnusComparisons": magnus_comparisons,
        "rowFlowCenters": row_flow_summaries,
        "rowFlowComparisons": row_flow_comparisons,
        "grassmannCenters": grassmann_summaries,
        "grassmannComparisons": grassmann_comparisons,
        "grassmannChartBallCertificate": chart_ball,
        "riccatiFlowEnclosure": riccati,
        "riccatiKrawczykCollocation": krawczyk,
        "coefficientBallEnclosure": coefficient_ball,
        "sampledCenterStabilityStatus": status(
            "sampled finite-difference endpoint center stability",
            sampled_stable,
            (
                "Stable only if mesh-refinement endpoint maps remain within a "
                "small constant factor.  They currently vary by many orders "
                "of magnitude."
            ),
        ),
        "exactSegmentCenterStabilityStatus": status(
            "exact-segment derivative endpoint center stability",
            exact_stable,
            (
                "Stable only if exact-derivative segment refinements remain "
                "within a small constant factor.  Current piecewise-constant "
                "propagation still varies by many orders of magnitude."
            ),
        ),
        "rankBallCenterReadyStatus": status(
            "rank-ball center ready",
            center_ready or cheb_stable or magnus_stable or row_flow_stable,
            (
                "A rank-ball proof needs a stable center before interval "
                "radii can be meaningful.  Current sampled and exact-segment "
                "centers are not stable enough; Chebyshev centers are tracked "
                "separately."
            ),
        ),
        "chebyshevCenterStabilityStatus": status(
            "Chebyshev endpoint center stability",
            cheb_stable,
            (
                "Closed only if all recorded Chebyshev refinement differences "
                "fall below the rank margin."
            ),
        ),
        "magnusCenterStabilityStatus": status(
            "Magnus endpoint center stability",
            magnus_stable,
            (
                "Closed only if all recorded Magnus refinement differences "
                "fall below the rank margin."
            ),
        ),
        "rowFlowMinorStabilityStatus": status(
            "adjoint row-flow minor stability",
            row_flow_stable,
            (
                "Closed only if the same transported-row 2x2 minor persists "
                "and changes by less than 25 percent across refinements."
            ),
        ),
        "grassmannProjectiveStabilityStatus": status(
            "normalized Grassmann projective stability",
            grassmann_stable,
            (
                "Closed at the numerical-center level if a persistent Pluecker "
                "chart stays away from zero and the tail projective comparisons "
                "are below tolerance.  This is not yet an interval enclosure "
                "of the exact endpoint flow."
            ),
        ),
        "grassmannChartBallStatus": status(
            "persistent Pluecker chart ball theorem",
            chart_ball_closed,
            (
                "Closed only as interval algebra: if the exact normalized "
                "Pluecker vector is enclosed in the stated projective ball, "
                "then the chart minor stays bounded away from zero."
            ),
        ),
        "riccatiFlowEquationStatus": status(
            "exact exterior/Riccati flow equation and tube",
            riccati_equations_closed,
            (
                "Closed when the exterior-square ODE, Riccati chart equation, "
                "and chart-tube noncollapse check are all recorded.  This "
                "still does not prove the exact endpoint BVP flow lies in the "
                "tube."
            ),
        ),
        "riccatiKrawczykCollocationStatus": status(
            "Riccati-chart Krawczyk collocation reduction",
            krawczyk_closed,
            (
                "Closed as a finite interval-collocation reduction.  Exact "
                "coefficient and endpoint-row interval balls still need to be "
                "proved below the reported Krawczyk capacities."
            ),
        ),
        "coefficientBallEnclosureStatus": status(
            "coefficient and boundary-row ball below Krawczyk capacity",
            coefficient_ball_closed,
            (
                "Closed under the stated sample/refinement ball model when "
                "the scaled companion and active boundary-row radii are below "
                "the Krawczyk capacities."
            ),
        ),
        "conclusion": (
            "The desired endpoint-map rank proof cannot be closed from the "
            "current sampled-flow center.  Chebyshev collocation greatly "
            "improves the scale and gives rank-two centers, but the recorded "
            "5->7 and 7->9 differences remain above the rank margin.  The "
            "two-point Magnus step is worse, blowing up under 2->4 "
            "refinement.  Adjoint row-flow propagation identifies the same "
            "best minor across refinements, but the raw determinant is not "
            "stable.  The normalized exterior/Grassmann diagnostic removes "
            "the raw scale and finds a persistent affine chart.  The "
            "chart-ball algebra now gives the exact projective radius needed "
            "for noncollapse.  The exterior/Riccati equations have been "
            "constructed, the Krawczyk/collocation reduction gives explicit "
            "coefficient and endpoint-row radius budgets, and the finite "
            "sample-ball enclosure is below those budgets with large slack. "
            "The last numerical rigor item is replacing the sample-ball model "
            "with analytic interval quadrature tails for the confluent "
            "integrals, or returning to obstruction annihilation."
        ),
        "nextTarget": (
            "Replace the coefficient sample-ball model with formal interval "
            "quadrature/tail bounds for the confluent integrals."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Endpoint map center stability")
    print(f"  sampled Frobenius spread: {sampled_spread:.12e}" if sampled_spread else "  sampled Frobenius spread: n/a")
    print(f"  exact-segment Frobenius spread: {exact_spread:.12e}" if exact_spread else "  exact-segment Frobenius spread: n/a")
    if cheb_comparisons:
        worst_cheb = max(float(row["differenceOverMargin"]) for row in cheb_comparisons)
        print(f"  worst Chebyshev diff/margin: {worst_cheb:.12e}")
    if magnus_comparisons:
        worst_magnus = max(float(row["differenceOverMargin"]) for row in magnus_comparisons)
        print(f"  worst Magnus diff/margin: {worst_magnus:.12e}")
    if row_flow_comparisons:
        worst_row_minor = max(
            float(row["bestMinorRelativeDiff"])
            for row in row_flow_comparisons
            if row.get("bestMinorRelativeDiff") is not None
        )
        print(f"  worst row-flow best-minor rel diff: {worst_row_minor:.12e}")
    if grassmann_comparisons:
        worst_grassmann = max(
            float(row["projectiveChordalDistance"])
            for row in grassmann_comparisons
        )
        print(f"  worst Grassmann projective distance: {worst_grassmann:.12e}")
    if chart_ball:
        print(
            "  chart-ball certified radius: "
            f"{float(chart_ball.get('certifiedProjectiveRadius', 0)):.12e}"
        )
        print(
            "  chart-ball threshold capacity: "
            f"{float(chart_ball.get('radiusCapacityForThreshold', 0)):.12e}"
        )
    if riccati:
        print(
            "  Riccati Lipschitz bound: "
            f"{float(riccati.get('riccatiLipschitzInfBoundOnChartTube', 0)):.12e}"
        )
    if krawczyk:
        print(
            "  Krawczyk coeff radius capacity: "
            f"{float(krawczyk.get('coefficientRadiusCapacityScaledCompanion', 0)):.12e}"
        )
        print(
            "  Krawczyk boundary radius capacity: "
            f"{float(krawczyk.get('boundaryRadiusCapacityUnscaledRows', 0)):.12e}"
        )
    if coefficient_ball:
        print(
            "  coefficient ball radius/capacity: "
            f"{float(coefficient_ball.get('scaledCompanionRadiusOverCapacity', 0)):.12e}"
        )
        print(
            "  boundary ball radius/capacity: "
            f"{float(coefficient_ball.get('boundaryRowRadiusOverCapacity', 0)):.12e}"
        )
    print(f"  sampled center stable: {sampled_stable}")
    print(f"  exact-segment center stable: {exact_stable}")
    print(f"  Chebyshev center stable: {cheb_stable}")
    print(f"  Magnus center stable: {magnus_stable}")
    print(f"  row-flow minor stable: {row_flow_stable}")
    print(f"  Grassmann projective stable: {grassmann_stable}")
    print(f"  chart-ball algebra closed: {chart_ball_closed}")
    print(f"  Riccati equations/tube closed: {riccati_equations_closed}")
    print(f"  Krawczyk collocation reduction closed: {krawczyk_closed}")
    print(f"  coefficient ball enclosure closed: {coefficient_ball_closed}")
    print(f"  rank-ball center ready: {center_ready or cheb_stable or magnus_stable or row_flow_stable}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
