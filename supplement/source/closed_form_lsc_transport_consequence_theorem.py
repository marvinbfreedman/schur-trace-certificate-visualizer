#!/usr/bin/env python3
r"""Narrow closed-form LSC transport consequence theorem.

This wrapper exposes only the closed lower-envelope facts needed by the
augmented nonnegative closed-form limit theorem.
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
        "--closed-form-lsc-json",
        default="closed_form_lsc_transport_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="closed_form_lsc_transport_consequence_theorem.json",
    )
    args = parser.parse_args()

    lsc = load(args.closed_form_lsc_json)
    envelope_ok = bool(lsc.get("closedLowerEnvelopeIdentifiedClosed"))
    transport_ok = bool(lsc.get("closedFormTransportClosed") or lsc.get("lowerSemicontinuityClosed"))
    ok = envelope_ok and transport_ok

    data = {
        "theoremName": "closed-form LSC transport consequence theorem",
        "proofClass": "symbolic identity",
        "closedFormLscJson": args.closed_form_lsc_json,
        "statement": (
            "The closed-form lower-semicontinuity transport theorem supplies "
            "the closed lower-envelope identification and transported LSC "
            "form facts needed by the augmented nonnegative closed-form limit theorem."
        ),
        "statuses": {
            "closedLowerEnvelopeIdentificationStatus": status(
                "closed lower-envelope identification",
                envelope_ok,
                "Imported from the closed-form LSC transport theorem.",
            ),
            "closedFormTransportStatus": status(
                "closed-form LSC transport",
                transport_ok,
                "Imported from the closed-form LSC transport theorem.",
            ),
            "closedFormLscTransportConsequenceStatus": status(
                "closed-form LSC transport consequence",
                ok,
                (
                    "Together these facts identify the completed lower "
                    "envelope used by the nonnegative closed-form limit theorem."
                ),
            ),
        },
        "closedLowerEnvelopeIdentifiedClosed": envelope_ok,
        "closedLowerEnvelopeLscClosed": envelope_ok,
        "lowerSemicontinuityClosed": transport_ok,
        "closedFormTransportClosed": transport_ok,
        "closedFormLscTransportConsequenceClosed": ok,
        "proof": [
            "Import the closed-form LSC transport theorem.",
            "Extract the closed lower-envelope identification.",
            "Extract the transported closed lower-semicontinuity conclusion.",
            "Expose only these consequences to the augmented nonnegative closed-form limit theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Closed lower-envelope identification or LSC transport is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed-form LSC transport consequence theorem")
    print(f"  closed lower-envelope identification: {envelope_ok}")
    print(f"  closed-form LSC transport: {transport_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
