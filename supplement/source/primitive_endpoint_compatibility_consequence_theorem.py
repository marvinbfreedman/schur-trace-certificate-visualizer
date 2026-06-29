#!/usr/bin/env python3
r"""Slim primitive endpoint compatibility consequence theorem."""

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
    parser.add_argument("--boundary-zero-json", default="primitive_boundary_zero_consequence_theorem.json")
    parser.add_argument("--density-json", default="primitive_trace_density_consequence_theorem.json")
    parser.add_argument("--dq-json", default="green_feature_dq_vanishes_theorem.json")
    parser.add_argument("--json-out", default="primitive_endpoint_compatibility_consequence_theorem.json")
    args = parser.parse_args()

    boundary = load(args.boundary_zero_json)
    density = load(args.density_json)
    dq = load(args.dq_json)

    boundary_ok = bool(boundary.get("D_bdyZeroOnPrimitiveTransport"))
    density_ok = bool(density.get("primitiveTraceDenseInXR")) and bool(
        density.get("DqZeroTransfersToPrimitiveTraceImage")
    )
    dq_ok = bool(dq.get("DqVanishesOnXR"))
    ok = boundary_ok and density_ok and dq_ok

    data = {
        "theoremName": "primitive endpoint compatibility consequence theorem",
        "proofClass": "symbolic identity",
        "boundaryZeroJson": args.boundary_zero_json,
        "densityJson": args.density_json,
        "dqJson": args.dq_json,
        "statuses": {
            "boundaryZeroInputStatus": status(
                "primitive boundary correction is zero",
                boundary_ok,
                "The primitive boundary-zero consequence gives D_bdy=0.",
            ),
            "dqZeroInputStatus": status(
                "quotient repair vanishes on X_R",
                dq_ok,
                "The Green-feature theorem gives D_q=0 on the completed trace range.",
            ),
            "dqVanishesOnXRStatus": status(
                "quotient repair D_q vanishes on X_R",
                dq_ok,
                (
                    "This is the exact trace-side consequence consumed by the "
                    "canonical boundary repair comparison: D_q=0 on the "
                    "completed transported trace range X_R."
                ),
            ),
            "primitiveDensityInputStatus": status(
                "D_q zero transfers to primitive trace image",
                density_ok,
                "The primitive trace image is dense in X_R and D_q is bounded.",
            ),
            "primitiveEndpointCompatibilityStatus": status(
                "primitive endpoint compatibility",
                ok,
                "The primitive boundary correction and quotient repair both vanish on primitive traces.",
            ),
        },
        "primitiveEndpointCompatibilityConsequenceClosed": ok,
        "primitiveEndpointCompatibilityClosed": ok,
        "proof": [
            "D_bdy=0 on primitive transport.",
            "D_q=0 on X_R by the Green-feature contraction.",
            "Primitive trace density transfers D_q=0 to primitive traces.",
            "Therefore endpoint/trace correction is annihilated for the original primitive lift.",
        ],
        "remainingAnalyticGap": None if ok else "Boundary-zero, D_q-zero, or density input is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive endpoint compatibility consequence theorem")
    print(f"  primitive endpoint compatibility: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
