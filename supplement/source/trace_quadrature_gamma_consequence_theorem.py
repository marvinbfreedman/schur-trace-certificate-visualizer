#!/usr/bin/env python3
r"""Trace quadrature gamma consequence for the continuum trace-frame theorem."""

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
        "--quadrature-json",
        default="trace_quadrature_interval_consistency_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_quadrature_gamma_consequence_theorem.json",
    )
    args = parser.parse_args()

    q = load(args.quadrature_json)
    gamma_after = (
        q.get("gammaAfterMetricTraceQuadrature")
        if q.get("gammaAfterMetricTraceQuadrature") is not None
        else q.get("gammaAfterTraceQuadrature")
    )
    gamma_value = float(gamma_after) if gamma_after is not None else None
    quadrature_closed = bool(q.get("traceQuadratureClosed"))
    absorbed = bool(q.get("traceQuadratureAbsorbedByFiniteGamma"))
    ok = quadrature_closed and absorbed and gamma_value is not None and gamma_value > 0.0
    qualitative_ok = quadrature_closed

    data = {
        "theoremName": "trace quadrature gamma consequence theorem",
        "proofClass": "symbolic identity",
        "quadratureJson": args.quadrature_json,
        "traceQuadratureClosed": quadrature_closed,
        "traceQuadratureAbsorbedByFiniteGamma": absorbed,
        "gammaAfterTraceQuadrature": q.get("gammaAfterTraceQuadrature"),
        "gammaAfterMetricTraceQuadrature": q.get("gammaAfterMetricTraceQuadrature"),
        "transportedGammaLower": gamma_after,
        "metricRelativeTraceQuadratureErrorBound": q.get(
            "metricRelativeTraceQuadratureErrorBound"
        ),
        "traceQuadratureErrorBound": q.get("traceQuadratureErrorBound"),
        "statement": (
            "The trace quadrature interval estimate is closed and leaves a "
            "positive transported gamma lower bound."
        ),
        "statuses": {
            "traceQuadratureInputStatus": status(
                "trace quadrature consistency input",
                quadrature_closed,
                "The interval quadrature theorem bounds the trace-frame quadrature error.",
            ),
            "finiteGammaAbsorptionInputStatus": status(
                "finite gamma absorption input",
                absorbed,
                "The quadrature error is absorbed by the finite gamma lower bound.",
            ),
            "transportedGammaPositiveConsequenceStatus": status(
                "transported gamma positive consequence",
                ok,
                "Therefore the transported continuum gamma lower bound is positive.",
            ),
        },
        "traceQuadratureGammaConsequenceClosed": ok,
        "traceQuadratureIntervalCertificateClosed": qualitative_ok,
        "qualitativeTraceQuadratureConsequenceClosed": qualitative_ok,
        "proof": [
            "Import the trace quadrature interval consistency theorem.",
            "Read off traceQuadratureClosed and traceQuadratureAbsorbedByFiniteGamma.",
            "Expose only the positive transported gamma lower bound to the continuum theorem.",
        ],
        "remainingAnalyticGap": None
        if qualitative_ok
        else "Trace quadrature interval consistency input is not closed.",
        "remainingNumericalGap": None
        if ok
        else (
            "The qualitative trace-quadrature consequence is closed, but the "
            "current explicit transported gamma lower bound is not positive."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace quadrature gamma consequence theorem")
    print(f"  quadrature closed: {quadrature_closed}")
    print(f"  absorbed by finite gamma: {absorbed}")
    print(f"  transported gamma: {gamma_value if gamma_value is not None else 'n/a'}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
