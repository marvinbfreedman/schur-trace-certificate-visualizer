#!/usr/bin/env python3
"""Package the undivided graph/flow Krawczyk result without overpromotion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": bool(closed),
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--graph-json", default="endpoint_undivided_dae_krawczyk_enclosure.json"
    )
    parser.add_argument(
        "--flow-json", default="endpoint_undivided_flow_endpoint_krawczyk.json"
    )
    parser.add_argument(
        "--crossing-json", default="endpoint_crossing_continuum_promotion_theorem.json"
    )
    parser.add_argument(
        "--json-out", default="endpoint_undivided_dae_flow_interval_theorem.json"
    )
    args = parser.parse_args()

    graph = load(args.graph_json)
    flow = load(args.flow_json)
    crossing = load(args.crossing_json)
    graph_closed = bool(
        graph.get("finitePolynomialDaeKrawczykStatus", {}).get("closed")
    )
    flow_closed = bool(
        flow.get("finiteInducedFlowKrawczykStatus", {}).get("closed")
    )
    finite_minor_closed = bool(
        flow.get("finiteRestrictedEndpointMinorStatus", {}).get("closed")
    )
    coefficient_cover_closed = bool(
        graph.get("exactEigenrowCoefficientCoverStatus", {}).get("closed")
    )
    exact_graph_closed = bool(
        graph.get("validatedExactUndividedDaeFlowStatus", {}).get("closed")
    )
    exact_endpoint_closed = bool(
        flow.get("exactContinuumEndpointFlowStatus", {}).get("closed")
    )
    exact_crossing_closed = bool(
        crossing.get("exactCrossingContinuumPromotionClosed")
    )
    finite_closed = graph_closed and flow_closed and finite_minor_closed
    exact_closed = (
        finite_closed
        and exact_crossing_closed
        and coefficient_cover_closed
        and exact_graph_closed
        and exact_endpoint_closed
    )

    data = {
        "theoremName": "undivided DAE graph and endpoint-flow collocation interval theorem",
        "proofClass": "Arb interval collocation/Krawczyk theorem",
        "proofArtifact": "endpoint_undivided_dae_flow_interval_theorem.md",
        "dependencies": [args.graph_json, args.flow_json, args.crossing_json],
        "domain": graph.get("domain"),
        "equations": {
            "graph": "a8(q' + qB + (qb)q) + c7 q + d = 0",
            "inducedFlow": "G' = (B+bq)G",
            "endpointMap": "M = C_R T(q_R)G_R - C_L T(q_L)G_L",
        },
        "finitePolynomialUndividedFlowStatus": status(
            "undivided graph and induced flow in the named polynomial collocation model",
            graph_closed and flow_closed,
            (
                "Every graph and induced-flow collocation box has an outward-rounded "
                "Arb Krawczyk self-map, including the chart containing the zero of a8."
            ),
        ),
        "finitePolynomialEndpointRankStatus": status(
            "restricted active endpoint map has full row rank in the polynomial collocation model",
            finite_closed,
            (
                "The interval endpoint matrix has a 2 by 2 minor whose Arb ball "
                "excludes zero."
            ),
        ),
        "exactEigenrowCoefficientCoverStatus": graph.get(
            "exactEigenrowCoefficientCoverStatus", {}
        ),
        "exactCrossingContinuumPromotionStatus": crossing.get(
            "exactCrossingContinuumPromotionStatus", {}
        ),
        "exactContinuumPromotionStatus": status(
            "exact Riemann-kernel undivided flow and endpoint rank",
            exact_closed,
            (
                "The finite polynomial theorem and the exact singular-crossing promotion "
                "are closed. Full endpoint promotion still requires exact regular-piece "
                "propagation and certification of the restricted active endpoint minor."
            ),
        ),
        "bestFiniteEndpointMinor": flow.get("bestIntervalMinor"),
        "finitePolynomialCollocationTheoremClosed": finite_closed,
        "finitePolynomialTheoremClosed": finite_closed,
        "exactCrossingContinuumPromotionClosed": exact_crossing_closed,
        "exactContinuumTheoremClosed": exact_closed,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"finite polynomial theorem closed: {finite_closed}")
    print(f"exact continuum theorem closed: {exact_closed}")
    print(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
