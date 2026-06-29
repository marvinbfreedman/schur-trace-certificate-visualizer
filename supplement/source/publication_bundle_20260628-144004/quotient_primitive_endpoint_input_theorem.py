#!/usr/bin/env python3
r"""Primitive endpoint input for quotient-to-original lift."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primitive-endpoint-json",
        default="primitive_endpoint_compatibility_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="quotient_primitive_endpoint_input_theorem.json")
    args = parser.parse_args()

    primitive = load(args.primitive_endpoint_json)
    primitive_ok = bool(primitive.get("primitiveEndpointCompatibilityClosed"))
    data = {
        "theoremName": "quotient primitive endpoint input theorem",
        "proofClass": "symbolic identity",
        "primitiveEndpointJson": args.primitive_endpoint_json,
        "statuses": {
            "primitiveEndpointInputStatus": status(
                "primitive endpoint compatibility input",
                primitive_ok,
                "The primitive endpoint consequence supplies the endpoint compatibility condition needed by the quotient lift.",
            ),
        },
        "quotientPrimitiveEndpointInputClosed": primitive_ok,
        "primitiveEndpointCompatibilityClosed": primitive_ok,
        "proof": [
            "Import the primitive endpoint compatibility consequence.",
            "Expose it under the quotient-lift input name.",
        ],
        "remainingAnalyticGap": None if primitive_ok else "Primitive endpoint compatibility consequence is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Quotient primitive endpoint input theorem")
    print(f"  primitive endpoint input: {primitive_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
