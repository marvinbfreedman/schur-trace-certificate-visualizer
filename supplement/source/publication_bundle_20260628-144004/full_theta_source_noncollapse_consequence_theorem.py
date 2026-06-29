#!/usr/bin/env python3
r"""Terminal consequence for full-theta source noncollapse.

The high-block source spectral-gap theorem only needs the closed full-theta
source noncollapse flag, the Riesz projector stability status, and the positive
gap/margin numbers.  This wrapper exports exactly that contract.
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
    if not item and isinstance(data.get("statuses"), dict):
        item = data["statuses"].get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def positive(data: dict[str, Any], key: str) -> bool:
    try:
        return float(data[key]) > 0.0
    except Exception:
        return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-json",
        default="full_theta_source_noncollapse_interval_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="full_theta_source_noncollapse_consequence_theorem.json",
    )
    args = parser.parse_args()

    source = load(args.source_json)
    noncollapse_ok = bool(
        source.get("fullThetaSourceNoncollapseIntervalTheoremClosed")
        or source.get("fullPhiContinuumSourceNoncollapsePasses")
        or nested_closed(source, "fullThetaSourceNoncollapseStatus")
    )
    projector_ok = nested_closed(source, "rieszProjectorStabilityStatus")
    gap_ok = positive(source, "spectralGapToComplementS8")
    lower_ok = positive(source, "lowerBoundAfterContinuumAndTail")
    closed = noncollapse_ok and projector_ok and gap_ok and lower_ok

    data: dict[str, Any] = {
        "theoremName": "full-theta source noncollapse terminal consequence theorem",
        "proofClass": "symbolic identity",
        "omega": source.get("omega"),
        "L": source.get("L"),
        "basis": source.get("basis"),
        "activeDimension": source.get("activeDimension"),
        "activeEigenvaluesS8": source.get("activeEigenvaluesS8", []),
        "spectralGapToComplementS8": source.get("spectralGapToComplementS8"),
        "lowerBoundAfterContinuumAndTail": source.get(
            "lowerBoundAfterContinuumAndTail"
        ),
        "totalContinuumErrorBound": source.get("totalContinuumErrorBound"),
        "totalProjectorAlpha": source.get("totalProjectorAlpha"),
        "fullThetaSourceNoncollapseIntervalTheoremClosed": closed,
        "fullPhiContinuumSourceNoncollapsePasses": closed,
        "rieszProjectorStabilityStatus": status(
            "active Riesz projector stability",
            projector_ok,
            "The full source noncollapse theorem closes the Riesz projector stability status.",
        ),
        "fullThetaSourceNoncollapseStatus": status(
            "full-Phi active source noncollapse",
            noncollapse_ok,
            "The full source noncollapse theorem proves the full-Phi active rank survives.",
        ),
        "positiveGapStatus": status(
            "positive source spectral gap and post-tail lower bound",
            gap_ok and lower_ok,
            "The imported source theorem reports positive spectral gap and post-tail lower response margin.",
        ),
        "fullThetaSourceNoncollapseConsequenceClosed": closed,
        "proof": [
            "Import the full-theta source noncollapse interval theorem.",
            "Read off the Riesz projector stability status and positive gap/margin scalars.",
            "Export only the source spectral-gap contract needed by the high-block source spectral-gap theorem.",
        ],
        "remainingAnalyticGap": None
        if closed
        else "Full-theta source noncollapse, Riesz stability, or positive gap is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Full-theta source noncollapse consequence theorem")
    print(f"  noncollapse: {noncollapse_ok}")
    print(f"  Riesz stability: {projector_ok}")
    print(f"  positive gap: {gap_ok}")
    print(f"  positive lower margin: {lower_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
