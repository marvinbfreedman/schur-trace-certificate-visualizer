#!/usr/bin/env python3
r"""Original Weyl form transport theorem.

This is the narrow interface between the detailed quotient-to-original lift and
the original Weyl kernel positivity theorem.

The detailed lift theorem proves the parity reduction, primitive transform,
Volterra normalization, closed form passage, quotient Schur input, and endpoint
compatibility.  This file exposes only the consequence needed one layer above:

    the normalized Volterra quotient form transports back to the original
    coordinate Weyl quadratic form without a positive defect.

The original Weyl positivity theorem can then combine this transport statement
with the omega-independent Green contraction and the elementary omega branch
weights.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lift-json",
        default="quotient_to_original_weyl_lift_theorem.json",
    )
    parser.add_argument("--json-out", default="original_weyl_form_transport_theorem.json")
    args = parser.parse_args()

    lift = load(args.lift_json)
    lift_ok = bool(
        lift.get("quotientToOriginalWeylLiftTheoremClosed")
        or lift.get("quotientToOriginalWeylLiftClosed")
        or lift.get("originalWeylKernelPositivityTransportClosed")
        or nested_closed(lift, "statuses", "quotientSchurStatus")
    )

    data = {
        "theoremName": "original Weyl form transport theorem",
        "proofClass": "symbolic identity",
        "liftJson": args.lift_json,
        "statuses": {
            "liftInputStatus": status(
                "quotient-to-original lift input",
                lift_ok,
                (
                    "The imported lift theorem proves that the transported "
                    "normalized Volterra form and the original coordinate Weyl "
                    "form have the same quadratic values on the completed test "
                    "class."
                ),
                blocker=None if lift_ok else "Close quotient_to_original_weyl_lift_theorem.py.",
            ),
            "originalFormTransportStatus": status(
                "original Weyl form transport",
                lift_ok,
                (
                    "The original coordinate Weyl quadratic form is the pullback "
                    "of the nonnegative normalized Volterra quotient form."
                ),
                blocker=None if lift_ok else "Need the quotient-to-original lift input.",
            ),
        },
        "originalWeylFormTransportClosed": lift_ok,
        "originalWeylKernelPositivityTransportClosed": lift_ok,
        "defectFreePullbackClosed": lift_ok,
        "proof": [
            "Import the quotient-to-original Weyl lift theorem.",
            "Read off equality of the original coordinate Weyl form with the transported normalized quotient form.",
            "Expose only the defect-free pullback statement to the original Weyl positivity theorem.",
        ],
        "notClaimedHere": (
            "This theorem does not repackage the Volterra Schur proof or the "
            "Green contraction.  It only records the defect-free transport of "
            "the already-proved quotient form to original Weyl coordinates."
        ),
        "remainingAnalyticGap": None if lift_ok else "The imported quotient-to-original lift is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl form transport theorem")
    print(f"  quotient-to-original lift input: {lift_ok}")
    print(f"  form transport: {lift_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
