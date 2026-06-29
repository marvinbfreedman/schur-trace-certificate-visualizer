#!/usr/bin/env python3
r"""Diagnostic consequence of the continuum trace-to-source Green kernel ledger.

This wrapper exposes only the supporting trace-to-source diagnostics consumed
by canonical boundary repair comparison.  It does not claim the full continuum
trace-to-source Green theorem is closed.
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


def closed_status(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--trace-green-json",
        default="continuum_trace_to_source_green_kernel.json",
    )
    parser.add_argument(
        "--json-out",
        default="continuum_trace_to_source_green_diagnostic_consequence_theorem.json",
    )
    args = parser.parse_args()

    trace_green = load(args.trace_green_json)
    finite_range_ok = closed_status(trace_green, "finiteActiveRangeIdentityStatus")
    sampled_density_ok = closed_status(trace_green, "sampledGreenDensityStatus")
    density_stability_ok = closed_status(trace_green, "traceGridDensityStabilityStatus")
    lagrange_ok = closed_status(trace_green, "lagrangeAdjointIdentityStatus")
    diagnostic_ok = finite_range_ok and sampled_density_ok and density_stability_ok and lagrange_ok

    data = {
        "theoremName": "continuum trace-to-source Green diagnostic consequence theorem",
        "proofClass": "symbolic identity",
        "traceGreenJson": args.trace_green_json,
        "statement": (
            "The continuum trace-to-source Green kernel ledger supplies finite "
            "active-range, sampled-density, density-stability, and Lagrange "
            "identity diagnostics used only as imported signals by the canonical "
            "boundary comparison.  The full continuum trace-to-source theorem "
            "remains open below this diagnostic consequence."
        ),
        "statuses": {
            "finiteActiveRangeIdentityStatus": status(
                "finite active range identity",
                finite_range_ok,
                "Imported from the continuum trace-to-source Green kernel ledger.",
            ),
            "sampledGreenDensityStatus": status(
                "sampled Green density",
                sampled_density_ok,
                "Imported from the continuum trace-to-source Green kernel ledger.",
            ),
            "traceGridDensityStabilityStatus": status(
                "trace-grid density stability",
                density_stability_ok,
                "Imported from the continuum trace-to-source Green kernel ledger.",
            ),
            "lagrangeAdjointIdentityStatus": status(
                "Lagrange adjoint identity",
                lagrange_ok,
                "Imported from the continuum trace-to-source Green kernel ledger.",
            ),
            "traceToSourceDiagnosticConsequenceStatus": status(
                "trace-to-source diagnostic consequence",
                diagnostic_ok,
                (
                    "The finite/sampled diagnostic signals are available for "
                    "audit display without importing the open continuum Green "
                    "theorem as a closed dependency."
                ),
            ),
        },
        "finiteActiveRangeIdentityClosed": finite_range_ok,
        "sampledGreenDensityClosed": sampled_density_ok,
        "traceGridDensityStabilityClosed": density_stability_ok,
        "lagrangeAdjointIdentityClosed": lagrange_ok,
        "traceToSourceDiagnosticConsequenceClosed": diagnostic_ok,
        "continuumTraceToSourceClosed": bool(trace_green.get("continuumTraceToSourceClosed")),
        "remainingAnalyticGap": None
        if diagnostic_ok
        else "One trace-to-source diagnostic input is open.",
        "proof": [
            "Import the continuum trace-to-source Green kernel ledger.",
            "Extract only closed finite/sampled diagnostic signals used for importedSignals.",
            "Do not promote the open continuum trace-to-source theorem.",
        ],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum trace-to-source Green diagnostic consequence theorem")
    print(f"  finite active range identity: {finite_range_ok}")
    print(f"  sampled Green density: {sampled_density_ok}")
    print(f"  trace-grid density stability: {density_stability_ok}")
    print(f"  Lagrange identity: {lagrange_ok}")
    print(f"  consequence: {diagnostic_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
