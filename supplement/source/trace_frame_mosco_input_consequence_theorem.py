#!/usr/bin/env python3
r"""Narrow trace-frame Mosco input consequence theorem.

This wrapper exposes only the model-specific trace-frame facts needed by the
high-block Mosco projection input theorem:

    continuum trace-frame lower bound,
    trace quadrature consistency,
    bounded sampled-trace correction right inverse.
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


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    if not item and isinstance(data.get("statuses"), dict):
        item = data["statuses"].get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--trace-frame-json",
        default="continuum_trace_frame_lower_bound_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_frame_mosco_input_consequence_theorem.json",
    )
    args = parser.parse_args()

    trace = load(args.trace_frame_json)
    lower_ok = closed(trace, "continuumTraceFrameLowerBoundStatus")
    quadrature_ok = closed(trace, "traceQuadratureConsistencyStatus")
    correction_ok = closed(trace, "boundedSampledTraceCorrectionStatus")
    ok = lower_ok and quadrature_ok and correction_ok

    data = {
        "theoremName": "trace-frame Mosco input consequence theorem",
        "proofClass": "symbolic identity",
        "traceFrameJson": args.trace_frame_json,
        "explicitContinuumGammaLower": trace.get("explicitContinuumGammaLower"),
        "activeDimension": trace.get("activeDimension"),
        "statement": (
            "The continuum trace-frame lower-bound theorem supplies the "
            "model-specific trace-frame inputs needed for high-block Mosco "
            "projection convergence: a positive continuum trace-frame lower "
            "bound, quadrature consistency, and a bounded sampled correction "
            "right inverse."
        ),
        "statuses": {
            "continuumTraceFrameLowerBoundStatus": status(
                "continuum trace-frame lower bound",
                lower_ok,
                "Imported from the continuum trace-frame theorem.",
            ),
            "traceQuadratureConsistencyStatus": status(
                "trace quadrature consistency",
                quadrature_ok,
                "Imported from the continuum trace-frame theorem.",
            ),
            "boundedSampledTraceCorrectionStatus": status(
                "bounded sampled-trace correction",
                correction_ok,
                "Imported from the continuum trace-frame theorem.",
            ),
            "traceFrameMoscoInputConsequenceStatus": status(
                "trace-frame Mosco input consequence",
                ok,
                (
                    "Together these three facts give the model-specific "
                    "high-block recovery and closed-trace compatibility input."
                ),
            ),
        },
        "traceFrameMoscoInputClosed": ok,
        "continuumTraceFrameLowerBoundClosed": lower_ok,
        "traceQuadratureConsistencyClosed": quadrature_ok,
        "boundedSampledTraceCorrectionClosed": correction_ok,
        "proof": [
            "Import the continuum trace-frame lower bound.",
            "Import trace quadrature consistency on the active block.",
            "Import the bounded sampled-trace correction right inverse.",
            "Expose only these three inputs to the high-block Mosco projection theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Trace-frame lower bound, quadrature consistency, or sampled correction is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace-frame Mosco input consequence theorem")
    print(f"  continuum trace-frame lower bound: {lower_ok}")
    print(f"  trace quadrature consistency: {quadrature_ok}")
    print(f"  bounded sampled correction: {correction_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
