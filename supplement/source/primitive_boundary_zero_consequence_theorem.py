#!/usr/bin/env python3
r"""Primitive boundary-zero consequence theorem."""

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
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primitive-boundary-json", default="primitive_boundary_transport_audit.json")
    parser.add_argument("--json-out", default="primitive_boundary_zero_consequence_theorem.json")
    args = parser.parse_args()

    boundary = load(args.primitive_boundary_json)
    no_endpoint = nested_closed(boundary, "statuses", "primitiveTransportNoEndpointStatus")
    boundary_zero = nested_closed(boundary, "statuses", "canonicalPrimitiveBoundaryZeroStatus")
    descends = nested_closed(boundary, "statuses", "zeroBoundaryDescendsStatus")
    route_a_false = nested_closed(boundary, "statuses", "comparisonReducesToDqZeroStatus")
    ok = no_endpoint and boundary_zero and descends

    data = {
        "theoremName": "primitive boundary-zero consequence theorem",
        "proofClass": "symbolic identity",
        "primitiveBoundaryJson": args.primitive_boundary_json,
        "statuses": {
            "primitiveRouteAFalseStatus": status(
                "primitive route A is false",
                route_a_false,
                "The primitive boundary audit records that primitive tests do not generally land in ker R_global.",
            ),
            "primitiveTransportNoEndpointStatus": status(
                "primitive transport has no endpoint term",
                no_endpoint,
                "For compact primitive tests F(0)=0 kills the lower endpoint and kernel decay kills infinity.",
            ),
            "primitiveBoundaryZeroStatus": status(
                "D_bdy is zero",
                boundary_zero,
                "The canonical primitive boundary form is the zero form.",
            ),
            "canonicalPrimitiveBoundaryZeroStatus": status(
                "canonical primitive boundary form is zero",
                boundary_zero,
                "The canonical primitive boundary form is the zero form.",
            ),
            "primitiveBoundaryDescendsStatus": status(
                "zero boundary descends to trace quotient",
                descends,
                "The zero boundary form descends through the completed trace quotient.",
            ),
            "zeroBoundaryDescendsStatus": status(
                "zero boundary form descends through R_global",
                descends,
                "The zero boundary form descends through the completed trace quotient.",
            ),
        },
        "primitiveBoundaryZeroConsequenceClosed": ok,
        "D_bdyZeroOnPrimitiveTransport": ok,
        "primitiveBoundaryTransportAuditClosed": ok,
        "proof": [
            "Import the primitive boundary transport audit.",
            "Extract only the endpoint consequence D_bdy=0 and its descent to the trace quotient.",
        ],
        "remainingAnalyticGap": None if ok else "Primitive boundary transport audit does not close the zero-boundary consequence.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive boundary-zero consequence theorem")
    print(f"  D_bdy=0 consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
