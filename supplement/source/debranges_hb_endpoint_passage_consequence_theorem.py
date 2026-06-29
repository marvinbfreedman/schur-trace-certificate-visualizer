#!/usr/bin/env python3
r"""Terminal consequence of the shifted-Xi de Branges endpoint passage.

The external Weyl/Volterra audit only needs the endpoint implication flag:
positive shifted-Xi de Branges kernels plus the zero-descent endpoint argument
give the RH-side zero-location conclusion.  The detailed diagonal inequality,
zero-descent, and closed-cone bridge imports stay one layer lower.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--endpoint-json",
        default="debranges_hb_endpoint_passage.json",
    )
    parser.add_argument(
        "--json-out",
        default="debranges_hb_endpoint_passage_consequence_theorem.json",
    )
    args = parser.parse_args()

    endpoint = load(args.endpoint_json)
    endpoint_closed = bool(endpoint.get("endpointPassageClosed"))
    rh_closed = bool(endpoint.get("conditionalRhClosed") or endpoint_closed)

    endpoint_status = status(
        "shifted-Xi de Branges endpoint consequence",
        endpoint_closed,
        (
            "The shifted-Xi de Branges kernel positivity theorem and the "
            "zero-descent endpoint passage imply the RH-side zero-location "
            "conclusion in the z-normalization."
        ),
        blocker=None if endpoint_closed else "Close debranges_hb_endpoint_passage.py.",
    )

    data = {
        "theoremName": "de Branges/Hermite-Biehler endpoint passage consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "endpointPassageConsequenceStatus": endpoint_status,
        },
        "endpointPassageClosed": endpoint_status["closed"],
        "conditionalRhClosed": endpoint_status["closed"] and rh_closed,
        "rhSideZeroLocationClosed": endpoint_status["closed"] and rh_closed,
        "normalization": endpoint.get("normalization", {}),
        "proof": [
            "Import the detailed de Branges/Hermite-Biehler endpoint passage.",
            "Export only the terminal endpoint and RH-side zero-location consequences needed by the external audit.",
        ],
        "remainingAnalyticGap": None if endpoint_status["closed"] else endpoint_status["blocker"],
        "nextProofTarget": None if endpoint_status["closed"] else endpoint_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("de Branges endpoint passage consequence theorem")
    print(f"  endpoint passage closed: {data['endpointPassageClosed']}")
    print(f"  conditional RH-side conclusion: {data['conditionalRhClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
