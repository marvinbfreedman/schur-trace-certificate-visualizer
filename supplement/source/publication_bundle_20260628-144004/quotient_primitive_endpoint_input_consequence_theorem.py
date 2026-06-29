#!/usr/bin/env python3
r"""Slim consequence for the quotient primitive endpoint input.

The publication quotient-to-original lift theorem only needs the closed input
flag saying primitive endpoint compatibility is available for the quotient
lift.  This wrapper exposes that single flag and leaves the primitive endpoint
compatibility chain below ``quotient_primitive_endpoint_input_theorem``.
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
        "--endpoint-input-json",
        default="quotient_primitive_endpoint_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="quotient_primitive_endpoint_input_consequence_theorem.json",
    )
    args = parser.parse_args()

    endpoint = load(args.endpoint_input_json)
    ok = bool(
        endpoint.get("quotientPrimitiveEndpointInputClosed")
        or endpoint.get("primitiveEndpointCompatibilityClosed")
    )

    data = {
        "theoremName": "quotient primitive endpoint input consequence theorem",
        "proofClass": "symbolic identity",
        "endpointInputSource": "quotient primitive endpoint input theorem",
        "statuses": {
            "primitiveEndpointInputStatus": status(
                "primitive endpoint compatibility input for quotient lift",
                ok,
                "The imported quotient primitive endpoint input theorem closes the endpoint compatibility input used by the lift.",
            )
        },
        "quotientPrimitiveEndpointInputClosed": ok,
        "primitiveEndpointCompatibilityClosed": ok,
        "proof": [
            "Import the quotient primitive endpoint input theorem.",
            "Expose only the closed endpoint-input flag consumed by the publication lift theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Quotient primitive endpoint input theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Quotient primitive endpoint input consequence theorem")
    print(f"  endpoint input: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
