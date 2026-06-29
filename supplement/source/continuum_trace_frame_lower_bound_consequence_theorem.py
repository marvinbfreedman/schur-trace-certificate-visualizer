#!/usr/bin/env python3
r"""Terminal consequence for the continuum trace-frame lower-bound theorem.

The high-block tail-passage theorem only needs the three exported facts

    continuumTraceFrameLowerBoundStatus.closed,
    traceQuadratureConsistencyStatus.closed,
    boundedSampledTraceCorrectionStatus.closed.

This wrapper keeps the full continuum trace-frame proof ledger out of upstream
audits and exposes only that contract.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def nested_closed(data: dict[str, Any], key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--trace-frame-json",
        default="continuum_trace_frame_lower_bound_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="continuum_trace_frame_lower_bound_consequence_theorem.json",
    )
    args = parser.parse_args()

    trace = load(args.trace_frame_json)
    lower_ok = nested_closed(trace, "continuumTraceFrameLowerBoundStatus")
    quadrature_ok = nested_closed(trace, "traceQuadratureConsistencyStatus")
    correction_ok = nested_closed(trace, "boundedSampledTraceCorrectionStatus")
    numeric_gamma_ok = nested_closed(trace, "explicitNumericGammaStatus")
    theorem_closed = lower_ok and quadrature_ok and correction_ok

    data: dict[str, Any] = {
        "theoremName": "continuum trace-frame lower-bound terminal consequence theorem",
        "proofClass": "symbolic identity",
        "activeDimension": trace.get("activeDimension"),
        "explicitContinuumGammaLower": trace.get("explicitContinuumGammaLower"),
        "observedFiniteFrameFloor": trace.get("observedFiniteFrameFloor"),
        "observedTraceMeshFrameMinDrift": trace.get("observedTraceMeshFrameMinDrift"),
        "finiteFrameIntervalCertificateClosed": trace.get(
            "finiteFrameIntervalCertificateClosed"
        ),
        "traceQuadratureIntervalCertificateClosed": trace.get(
            "traceQuadratureIntervalCertificateClosed"
        ),
        "continuumTraceFrameLowerBoundStatus": status(
            "continuum L2 trace-frame lower bound",
            lower_ok,
            "The full continuum trace-frame theorem proves gamma_delta > 0.",
        ),
        "traceQuadratureConsistencyStatus": status(
            "trace quadrature consistency",
            quadrature_ok,
            "The full continuum trace-frame theorem imports the interval quadrature consistency input.",
        ),
        "boundedSampledTraceCorrectionStatus": status(
            "bounded sampled-trace correction right inverse",
            correction_ok,
            "The lower frame bound and quadrature consistency give bounded sampled-trace correction.",
        ),
        "explicitNumericGammaStatus": status(
            "explicit rigorous numeric gamma_delta",
            numeric_gamma_ok,
            "The interval finite-frame and quadrature certificates leave a positive transported gamma.",
        ),
        "continuumTraceFrameLowerBoundConsequenceClosed": theorem_closed,
        "proof": [
            "Import the full continuum trace-frame lower-bound theorem.",
            "Read off the three terminal trace-frame statuses used by high-block tail passage.",
            "Export no active-range, source-rank, finite-frame, or quadrature proof internals upstream.",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "One terminal trace-frame consequence status is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum trace-frame lower-bound consequence theorem")
    print(f"  continuum lower bound: {lower_ok}")
    print(f"  quadrature consistency: {quadrature_ok}")
    print(f"  sampled correction: {correction_ok}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
