#!/usr/bin/env python3
r"""Slim consequence for the Volterra/Green feature-map identity.

The Green-feature ``D_q=0`` theorem needs only:

* the Volterra moment representation of the trace Schur form on Green
  minimizer traces;
* the signed plus/minus square completion.

It does not use the open honest-positive-square question from the full feature
map ledger.  This wrapper exposes only the closed identification facts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


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
        "--feature-json",
        default="trace_volterra_green_feature_map.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_volterra_green_feature_consequence_theorem.json",
    )
    args = parser.parse_args()

    feature = load(args.feature_json)
    feature_def = nested_closed(feature, "statuses", "greenMinimizerFeatureDefinitionStatus")
    moment = nested_closed(feature, "statuses", "volterraMomentRepresentationStatus")
    signed_square = nested_closed(feature, "statuses", "signedSquareCompletionStatus")
    ok = feature_def and moment and signed_square

    data = {
        "theoremName": "Volterra/Green feature-map consequence theorem",
        "proofClass": "symbolic identity",
        "featureSource": "trace Volterra Green feature-map ledger",
        "featureDefinitions": feature.get("featureDefinitions"),
        "exactVolterraIdentity": feature.get("exactVolterraIdentity"),
        "signedFeatureIdentity": feature.get("signedFeatureIdentity"),
        "signedSquareIdentity": feature.get("signedSquareIdentity"),
        "statuses": {
            "greenMinimizerFeatureDefinitionStatus": status(
                "Green minimizer Volterra feature definition",
                feature_def,
                "The full feature-map ledger defines M, N, and G_+/G_- on Euler-Lagrange Green minimizers.",
            ),
            "volterraMomentRepresentationStatus": status(
                "Volterra moment representation of D_trace",
                moment,
                "The trace Schur form equals the explicit Volterra moment form on Green-minimizer traces.",
            ),
            "signedSquareCompletionStatus": status(
                "signed Volterra square completion",
                signed_square,
                "The elementary pointwise identity rewrites MN+NM as plus-square minus minus-square.",
            ),
        },
        "volterraGreenFeatureConsequenceClosed": ok,
        "greenFeatureIdentificationClosed": ok,
        "notClaimedHere": (
            "This consequence does not claim the open honest positive "
            "Volterra/Green Gram square from the full feature-map ledger."
        ),
        "proof": [
            "Import the closed Volterra moment representation on Green minimizers.",
            "Import the closed signed square completion identity.",
            "Expose only the feature-identification facts used by the D_q vanishing theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Feature definition, Volterra moment identity, or signed square completion is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Volterra/Green feature-map consequence theorem")
    print(f"  feature identification: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
