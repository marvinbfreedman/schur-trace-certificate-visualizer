#!/usr/bin/env python3
r"""Terminal consequence for source-inactive continuum tail absorption.

The global bridge only needs the absorption scalars and pass/fail flag from
the continuum tail absorption certificate.  This wrapper exports those fields
without preserving the raw inactive-tail certificate dependency upstream.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


PASS_KEYS = [
    "continuumTailEstimateProved",
    "absorptionCertificatePasses",
    "epsilonDeltaUpper",
    "finiteLowMidSchurBudget",
    "absorptionSlack",
    "epsilonOverBudget",
    "sourceGrid",
    "basis",
    "globalTraceRatio",
    "activeTraceKernelSourceFrac",
    "fullTraceKernelSourceFrac",
    "inactiveTraceKernelSourceFrac",
    "continuumInactiveTopUpper",
    "finiteInactiveTop",
    "continuumErrorBound",
    "conditionalTheorem",
    "remainingAnalyticGap",
]


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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
        "--absorption-json",
        default="continuum_tail_absorption_certificate.json",
    )
    parser.add_argument(
        "--json-out",
        default="continuum_tail_absorption_consequence_theorem.json",
    )
    args = parser.parse_args()

    absorption = load(args.absorption_json)
    summary = {key: absorption.get(key) for key in PASS_KEYS}
    absorption_pass = bool(summary["absorptionCertificatePasses"])
    try:
        slack_positive = float(summary["absorptionSlack"]) > 0.0
    except Exception:
        slack_positive = False
    consequence_closed = absorption_pass and slack_positive

    data: dict[str, Any] = {
        "theoremName": "source-inactive continuum tail absorption consequence theorem",
        "proofClass": "symbolic identity",
        **summary,
        "tailAbsorptionConsequenceStatus": status(
            "source-inactive tail absorption consequence",
            consequence_closed,
            (
                "The absorption certificate reports positive slack between "
                "the inactive epsilon upper bound and the finite low/mid "
                "Schur budget."
            ),
        ),
        "tailAbsorptionConsequenceClosed": consequence_closed,
        "proof": [
            "Import the continuum tail absorption certificate.",
            "Read off the absorption pass flag and scalar slack.",
            "Export only the absorption summary used by the global bridge.",
        ],
        "omittedRawPayloads": [
            "raw inactive-tail certificate pointer",
            "source/trace row table",
            "full source-inactive tail certificate payload",
        ],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum tail absorption consequence theorem")
    print(f"  absorption pass: {absorption_pass}")
    print(f"  slack positive: {slack_positive}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
