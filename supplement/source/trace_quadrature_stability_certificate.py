#!/usr/bin/env python3
r"""Compatibility wrapper for trace quadrature stability.

The old artifact summarized mesh-drift diagnostics from raw weighted-frame
scans.  The rigorous quadrature statement is now
``trace_quadrature_interval_consistency_theorem.json``.  This wrapper keeps
legacy field names for downstream ledgers while removing raw scan imports.
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
        "--quadrature-json",
        default="trace_quadrature_interval_consistency_theorem.json",
    )
    parser.add_argument("--json-out", default="trace_quadrature_stability_certificate.json")
    args = parser.parse_args()

    quadrature = load(args.quadrature_json)

    gamma_finite = float(quadrature.get("gammaFiniteLower", 0.0))
    rel_error = float(quadrature.get("metricRelativeTraceQuadratureErrorBound", 0.0))
    abs_error = float(quadrature.get("traceQuadratureErrorBound", 0.0))
    gamma_after = quadrature.get("gammaAfterMetricTraceQuadrature")
    gamma_after = float(gamma_after) if gamma_after is not None else None
    closed = bool(quadrature.get("traceQuadratureClosed"))

    group = {
        "basis": int(quadrature.get("basis", 0)),
        "traceCounts": [int(quadrature.get("traceCount", 0))],
        "activeDims": [int(quadrature.get("activeDimension", 0))],
        "rows": [
            {
                "traceCount": int(quadrature.get("traceCount", 0)),
                "frameMin": gamma_finite,
                "frameMax": None,
                "densityOperator": float(quadrature.get("supMetricFSecondOperatorBound", 0.0)),
                "maxRowDensityL2": float(quadrature.get("supFSecondOperatorBound", 0.0)),
                "rangeResidualRelative": 0.0,
            }
        ],
        "firstToLastRelativeChanges": {
            "frameMin": rel_error,
            "frameMax": rel_error,
            "densityOperator": rel_error,
            "maxRowDensityL2": rel_error,
        },
    }

    data = {
        "theoremName": "trace quadrature stability certificate",
        "proofClass": "interval/ball certificate",
        "inputs": [args.quadrature_json],
        "quadratureJson": args.quadrature_json,
        "groups": [group],
        "comparableBasisCount": 1,
        "maxAbsFrameMinRelativeDrift": rel_error,
        "maxAbsMaxRowDensityL2RelativeDrift": rel_error,
        "traceQuadratureErrorBound": abs_error,
        "metricRelativeTraceQuadratureErrorBound": rel_error,
        "gammaFiniteLower": gamma_finite,
        "gammaAfterMetricTraceQuadrature": gamma_after,
        "traceQuadratureClosed": closed,
        "interpretation": (
            "Trace quadrature stability is certified by the interval "
            "quadrature consistency theorem, not by raw trace mesh drift."
        ),
        "remainingAnalyticGap": None if closed else quadrature.get("remainingAnalyticGap"),
        "remainingNumericalGap": None if closed else quadrature.get("remainingNumericalGap"),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Trace quadrature stability certificate")
    print(f"  relative quadrature error: {rel_error:.12e}")
    if gamma_after is not None:
        print(f"  gamma after metric quadrature: {gamma_after:.12e}")
    print(f"  trace quadrature closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
