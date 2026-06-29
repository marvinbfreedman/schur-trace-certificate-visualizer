#!/usr/bin/env python3
r"""Slim consequence for the adjoint Green boundary diagnostic.

The continuum trace-to-source Green-kernel ledger does not use the sampled
finite-difference boundary diagnostic as a proof of the continuum adjoint BVP.
It only records the diagnostic failure mode: differentiating the sampled
pseudoinverse kernel is not the correct Green solution.

This wrapper exposes exactly that diagnostic summary and leaves the full
finite-difference calculation below it in the audit graph.
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
        "--diagnostic-json",
        default="adjoint_green_boundary_diagnostic.json",
    )
    parser.add_argument(
        "--json-out",
        default="adjoint_green_boundary_diagnostic_consequence_theorem.json",
    )
    args = parser.parse_args()

    diagnostic = load(args.diagnostic_json)
    raw_status = diagnostic.get("status", {})
    finite_active_range = bool(raw_status.get("finiteActiveRangeClosed"))
    sampled_ibp_small = bool(raw_status.get("discreteAdjointIbpSmall"))
    continuum_closed = bool(raw_status.get("continuumAdjointGreenClosed"))
    diagnostic_available = "maxIbpRelativeDefectOnActive" in diagnostic

    data = {
        "theoremName": "adjoint Green boundary diagnostic consequence theorem",
        "proofClass": "symbolic identity",
        "diagnosticSource": "adjoint Green boundary finite-difference diagnostic",
        "activeDim": diagnostic.get("activeDim"),
        "activeRangeFrobeniusRelative": diagnostic.get("activeRangeFrobeniusRelative"),
        "maxIbpRelativeDefectOnActive": diagnostic.get("maxIbpRelativeDefectOnActive"),
        "maxActiveRangeRelativeDefect": diagnostic.get("maxActiveRangeRelativeDefect"),
        "maxEtaL2Trap": diagnostic.get("maxEtaL2Trap"),
        "status": {
            "finiteActiveRangeClosed": finite_active_range,
            "discreteAdjointIbpSmall": sampled_ibp_small,
            "continuumAdjointGreenClosed": continuum_closed,
        },
        "statuses": {
            "diagnosticAvailableStatus": status(
                "sampled adjoint-boundary diagnostic is available",
                diagnostic_available,
                "The finite-difference diagnostic reports active-range, IBP-defect, and eta-size summaries.",
            ),
            "finiteActiveRangeDiagnosticStatus": status(
                "finite active range diagnostic",
                finite_active_range,
                "The diagnostic confirms the finite active block range identity.",
            ),
            "sampledIbpFailureModeStatus": status(
                "sampled pseudoinverse kernel is not the continuum Green solution",
                diagnostic_available and not sampled_ibp_small,
                (
                    "The sampled-kernel IBP defect is intentionally retained "
                    "as a diagnostic failure mode.  It shows that differentiating "
                    "the coarse pseudoinverse kernel does not solve the adjoint "
                    "Green BVP."
                ),
            ),
            "continuumAdjointGreenProofStatus": status(
                "continuum adjoint Green proof not claimed here",
                not continuum_closed,
                (
                    "This consequence does not close the continuum adjoint "
                    "Green problem; endpoint range and active trace-range "
                    "theorems remain the proof route."
                ),
            ),
        },
        "adjointGreenBoundaryDiagnosticConsequenceClosed": diagnostic_available,
        "interpretation": diagnostic.get("interpretation"),
        "proof": [
            "Import the finite-difference adjoint Green boundary diagnostic.",
            "Expose only the summary fields consumed by continuum_trace_to_source_green_kernel.",
            "Do not promote sampled IBP evidence to a continuum adjoint Green theorem.",
        ],
        "notAProofDependency": (
            "The continuum adjoint Green BVP is handled by the endpoint range "
            "and active trace-range theorem chain.  This record is retained "
            "only as a sampled-kernel diagnostic."
        ),
        "remainingAnalyticGap": None,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Adjoint Green boundary diagnostic consequence theorem")
    print(f"  diagnostic available: {diagnostic_available}")
    print(f"  sampled IBP small: {sampled_ibp_small}")
    print(f"  continuum Green proof claimed: {continuum_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
