#!/usr/bin/env python3
r"""Compatibility wrapper for the trace-frame continuum passage theorem.

The older version summarized raw weighted-frame scans.  The rigorous inputs
now live in:

  * trace_frame_interval_lower_bound_certificate.json
  * trace_quadrature_interval_consistency_theorem.json
  * continuum_trace_frame_lower_bound_theorem.json

This wrapper preserves the legacy fields consumed by high-block ledgers while
removing direct dependencies on raw scan JSON files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-frame-json",
        default="trace_frame_interval_lower_bound_certificate.json",
    )
    parser.add_argument(
        "--quadrature-json",
        default="trace_quadrature_interval_consistency_theorem.json",
    )
    parser.add_argument(
        "--continuum-json",
        default="continuum_trace_frame_lower_bound_theorem.json",
    )
    parser.add_argument("--json-out", default="trace_frame_continuum_passage_certificate.json")
    args = parser.parse_args()

    finite = load(args.finite_frame_json)
    quadrature = load(args.quadrature_json)
    continuum = load(args.continuum_json)

    gamma_finite = float(finite.get("gammaFiniteLower", 0.0))
    gamma_delta = continuum.get("explicitContinuumGammaLower")
    gamma_delta = float(gamma_delta) if gamma_delta is not None else None
    residual = float(finite.get("rangeResidualRelative", 0.0))
    finite_positive = bool(finite.get("gammaFiniteLowerPositive", gamma_finite > 0.0))
    quadrature_closed = bool(quadrature.get("traceQuadratureClosed"))
    continuum_closed = bool(
        continuum.get("remainingAnalyticGap") is None
        and continuum.get("remainingNumericalGap") is None
    )

    row = {
        "basis": int(finite.get("basis", 0)),
        "traceCount": int(finite.get("traceCount", 0)),
        "activeDim": int(finite.get("activeDimension", 0)),
        "frameRank": int(finite.get("activeDimension", 0)),
        "frameMin": gamma_finite,
        "frameMax": float(max(finite.get("centerFrameEigenvaluesFloat", [gamma_finite]))),
        "densityOperator": float(quadrature.get("supMetricFSecondOperatorBound", 0.0)),
        "maxRowDensityL2": float(quadrature.get("supFSecondOperatorBound", 0.0)),
        "frameBoundOperator": float(quadrature.get("metricRelativeTraceQuadratureErrorBound", 0.0)),
        "rangeResidualRelative": residual,
        "sourceFile": args.finite_frame_json,
    }

    data = {
        "theoremName": "trace-frame continuum passage certificate",
        "proofClass": "interval/ball certificate",
        "inputs": [
            args.finite_frame_json,
            args.quadrature_json,
            args.continuum_json,
        ],
        "finiteFrameJson": args.finite_frame_json,
        "quadratureJson": args.quadrature_json,
        "continuumTraceFrameJson": args.continuum_json,
        "rowCount": 1,
        "basisRange": [row["basis"], row["basis"]],
        "traceCounts": [row["traceCount"]],
        "activeDims": [row["activeDim"]],
        "allFrameMinsPositive": finite_positive,
        "observedGammaFloor": gamma_finite,
        "certifiedContinuumGammaLower": gamma_delta,
        "observedGammaFloorRow": row,
        "maxRangeResidualRelative": residual,
        "rangeResidualWithinTolerance": residual <= 1e-40,
        "maxDensityOperator": row["densityOperator"],
        "maxRowDensityL2": row["maxRowDensityL2"],
        "rows": [row],
        "finiteFrameIntervalCertificateClosed": finite_positive,
        "traceQuadratureIntervalCertificateClosed": quadrature_closed,
        "continuumTraceFrameLowerBoundClosed": continuum_closed,
        "continuumPassageClosed": finite_positive and quadrature_closed and continuum_closed,
        "conditionalConclusion": (
            "The continuum trace-frame lower bound and bounded sampled-trace "
            "correction right inverse are supplied by the interval lower-bound "
            "certificate, the interval quadrature theorem, and the continuum "
            "trace-frame lower-bound theorem."
        ),
        "remainingAnalyticGap": None if continuum_closed else continuum.get("remainingAnalyticGap"),
        "remainingNumericalGap": None if continuum_closed else continuum.get("remainingNumericalGap"),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Trace-frame continuum passage certificate")
    print(f"  finite gamma lower: {gamma_finite:.12e}")
    if gamma_delta is not None:
        print(f"  continuum gamma lower: {gamma_delta:.12e}")
    print(f"  interval quadrature closed: {quadrature_closed}")
    print(f"  continuum passage closed: {data['continuumPassageClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
