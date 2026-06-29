#!/usr/bin/env python3
r"""Primitive trace-density consequence theorem."""

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
    parser.add_argument("--primitive-density-json", default="primitive_trace_image_density.json")
    parser.add_argument("--json-out", default="primitive_trace_density_consequence_theorem.json")
    args = parser.parse_args()

    density = load(args.primitive_density_json)
    compact_core = nested_closed(density, "statuses", "primitiveContainsCompactCoreStatus")
    dense = nested_closed(density, "statuses", "primitiveTraceDenseStatus")
    dq_equiv = nested_closed(density, "statuses", "dqZeroEquivalenceStatus")
    ok = compact_core and dense and dq_equiv

    data = {
        "theoremName": "primitive trace-density consequence theorem",
        "proofClass": "symbolic identity",
        "primitiveDensityJson": args.primitive_density_json,
        "statuses": {
            "primitiveContainsCompactCoreStatus": status(
                "primitive image contains compact smooth core",
                compact_core,
                "Compact smooth Volterra core elements arise as derivatives of compact primitive tests.",
            ),
            "primitiveTraceDenseStatus": status(
                "primitive trace image is dense in X_R",
                dense,
                "The primitive trace image is dense in the completed transported trace range.",
            ),
            "dqZeroEquivalenceStatus": status(
                "D_q zero transfers from X_R to primitive image",
                dq_equiv,
                "Since D_q is bounded, D_q=0 on X_R is equivalent to D_q=0 on the dense primitive trace image.",
            ),
        },
        "primitiveTraceDensityConsequenceClosed": ok,
        "primitiveTraceDenseInXR": ok,
        "primitiveTraceImageDenseInXR": ok,
        "DqZeroTransfersToPrimitiveTraceImage": ok,
        "dqZeroOnPrimitiveImageEquivalentToDqZero": dq_equiv,
        "proof": [
            "Import the primitive trace image density theorem.",
            "Expose only the density and bounded-form transfer needed by endpoint compatibility.",
        ],
        "remainingAnalyticGap": None if ok else "Primitive trace density or D_q transfer input is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive trace-density consequence theorem")
    print(f"  primitive trace density consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
