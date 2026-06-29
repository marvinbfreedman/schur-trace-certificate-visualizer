#!/usr/bin/env python3
r"""Narrow consequence for the shifted-Xi zero-descent endpoint theorem.

The Hermite--Biehler endpoint passage only needs the terminal zero-location
facts:

    zeroDescentEndpointClosed = true,
    endpointZeroLocationClosed = true,
    conditionalRhClosed = true.

This file reads the full zero-descent theorem and exports only those endpoint
consequences.  The diagonal inequality and real-entire normalization inputs
remain below the full theorem.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


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
        "--zero-descent-json",
        default="shifted_xi_zero_descent_endpoint_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="shifted_xi_zero_descent_endpoint_consequence_theorem.json",
    )
    args = parser.parse_args()

    zero = load(args.zero_descent_json)
    zero_ok = bool(zero.get("zeroDescentEndpointClosed"))
    endpoint_ok = bool(zero.get("endpointZeroLocationClosed"))
    rh_ok = bool(zero.get("conditionalRhClosed"))
    ok = zero_ok and endpoint_ok and rh_ok

    data = {
        "theoremName": "shifted Xi zero-descent endpoint consequence theorem",
        "proofClass": "symbolic identity",
        "source": "shifted Xi zero-descent endpoint theorem",
        "statuses": {
            "zeroDescentEndpointInputStatus": status(
                "zero-descent endpoint input",
                zero_ok,
                "The full zero-descent theorem closes the endpoint contradiction.",
            ),
            "endpointZeroLocationInputStatus": status(
                "endpoint zero-location input",
                endpoint_ok,
                "The full zero-descent theorem proves all Xi zeros are real in z-normalization.",
            ),
            "conditionalRhInputStatus": status(
                "conditional RH-side input",
                rh_ok,
                "The full zero-descent theorem records the RH-side zero-location consequence.",
            ),
            "zeroDescentEndpointConsequenceStatus": status(
                "zero-descent endpoint consequence",
                ok,
                "Only the terminal zero-location consequence is exported upstream.",
            ),
        },
        "zeroDescentEndpointConsequenceClosed": ok,
        "zeroDescentEndpointClosed": zero_ok,
        "endpointZeroLocationClosed": endpoint_ok,
        "conditionalRhClosed": rh_ok,
        "proof": [
            "Import the full shifted-Xi zero-descent endpoint theorem.",
            "Export only the endpoint contradiction and zero-location flags used upstream.",
        ],
        "notExportedHere": [
            "diagonal de Branges inequality proof",
            "shifted-Xi normalization proof",
            "real-entire conjugation proof details",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Shifted-Xi zero-descent endpoint theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi zero-descent endpoint consequence theorem")
    print(f"  endpoint consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
