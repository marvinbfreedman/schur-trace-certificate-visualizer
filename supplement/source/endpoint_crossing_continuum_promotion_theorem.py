#!/usr/bin/env python3
"""Package the exact between-node crossing promotion theorem."""

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
        "--cover-json",
        default="endpoint_eigenrow_between_node_crossing_cover_q40.json",
    )
    parser.add_argument(
        "--abstract-json",
        default="continuous_undivided_dae_residual_inclusion_theorem.json",
    )
    parser.add_argument(
        "--application-json",
        default="endpoint_crossing_continuous_residual_inclusion_certificate.json",
    )
    parser.add_argument(
        "--json-out", default="endpoint_crossing_continuum_promotion_theorem.json"
    )
    args = parser.parse_args()

    cover = load(args.cover_json)
    abstract = load(args.abstract_json)
    application = load(args.application_json)
    cover_closed = bool(cover.get("betweenNodeEigenrowDerivativeCoverClosed"))
    abstract_closed = bool(abstract.get("continuousResidualInclusionTheoremClosed"))
    application_closed = bool(
        application.get("continuousCrossingResidualInclusionClosed")
    )
    closed = cover_closed and abstract_closed and application_closed

    data = {
        "schemaVersion": 1,
        "theoremName": "exact between-node promotion of the undivided endpoint crossing graph",
        "proofClass": "analytic proof plus interval/ball certificate",
        "proofArtifact": "endpoint_crossing_continuum_promotion_theorem.md",
        "statement": (
            "On the certified crossing interval, the exact normalized negative "
            "eigenrow derivatives e^(q), 0<=q<=8, lie in the displayed overlapping "
            "Bernstein/Cauchy cells, and the exact undivided graph DAE has a unique "
            "bounded solution in the certified scaled sup-norm tube."
        ),
        "domain": {
            "parameterInterval": cover.get("domain", {}).get("parameterInterval"),
            "state": "bounded seven-plane Riccati graph row in Taylor-scaled sup norm",
            "scope": (
                "exact crossing interval only; regular-piece propagation and the "
                "restricted active endpoint minor remain separate obligations"
            ),
        },
        "dependencies": [args.cover_json, args.abstract_json, args.application_json],
        "betweenNodeDerivativeCoverStatus": cover.get(
            "betweenNodeEigenrowDerivativeCoverStatus"
        ),
        "continuousResidualInclusionStatus": application.get(
            "continuousCrossingResidualInclusionStatus"
        ),
        "radiiPolynomial": application.get("radiiPolynomial"),
        "scaledConstants": application.get("scaledConstants"),
        "exactCrossingContinuumPromotionClosed": closed,
        "exactCrossingContinuumPromotionStatus": status(
            "exact continuum promotion through the simple endpoint crossing",
            closed,
            (
                "The exact derivative cover, singular inverse theorem, and concrete "
                "radii polynomial all close on the crossing interval."
                if closed
                else "At least one exact derivative-cover or residual-inclusion input remains open."
            ),
        ),
        "fullEndpointPropagationClosed": False,
        "fullEndpointPropagationStatus": status(
            "global exact propagation to both endpoint boundary maps",
            False,
            "This theorem closes the singular crossing only; it does not certify the two regular pieces or the final active endpoint minor.",
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"exact crossing promotion closed: {closed}")
    print(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
