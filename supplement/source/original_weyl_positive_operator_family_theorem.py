#!/usr/bin/env python3
r"""KLM-ready original Weyl positive-operator family theorem.

This is the narrow interface between the original Weyl kernel-to-operator
identity and the hbar=1 KLM packaging theorem.  The detailed
Volterra/Green/quotient-to-original machinery stays lower in the audit graph.
This file records only the operator-level consequence needed one layer above:

    for every |omega|<1/2, Op^W(sigma_omega) is positive semidefinite
    in the fixed hbar=1 Weyl normalization.

The conversion from this Weyl operator statement to KLM quantum positive type
is still handled separately by ``uniform_omega_weyl_klm_bridge.py`` using
``klm_weyl_hbar1_equivalence_theorem.json``.
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
        "--kernel-to-operator-json",
        default="original_weyl_kernel_to_operator_identity_theorem.json",
    )
    parser.add_argument(
        "--original-weyl-kernel-json",
        default="",
    )
    parser.add_argument(
        "--json-out",
        default="original_weyl_positive_operator_family_theorem.json",
    )
    args = parser.parse_args()

    identity = load(args.kernel_to_operator_json)
    original = load(args.original_weyl_kernel_json) if args.original_weyl_kernel_json else {}
    original_closed = bool(
        identity.get("originalWeylKernelToOperatorIdentityClosed")
        or identity.get("originalWeylOperatorPositiveClosed")
        or
        original.get("originalWeylKernelPositivityClosed")
        or nested_closed(original, "statuses", "originalWeylPositivityStatus")
    )

    data = {
        "theoremName": "original Weyl positive-operator family theorem",
        "proofClass": "symbolic identity",
        "omegaRange": "|omega| < 1/2",
        "normalization": "hbar=1 Weyl convention",
        "kernelToOperatorJson": args.kernel_to_operator_json,
        "legacyInputs": {"originalWeylKernelJson": args.original_weyl_kernel_json or None},
        "statuses": {
            "kernelToOperatorIdentityInputStatus": status(
                "kernel-to-operator identity input",
                original_closed,
                (
                    "The imported identity theorem converts original Weyl "
                    "kernel positivity into positive Weyl operators."
                ),
                blocker=None
                if original_closed
                else "Close original_weyl_kernel_to_operator_identity_theorem.py.",
            ),
            "weylOperatorFamilyStatus": status(
                "positive Weyl operator family",
                original_closed,
                (
                    "In the same Weyl normalization, positivity of the "
                    "coordinate kernel quadratic form is the statement that "
                    "Op^W(sigma_omega) is positive semidefinite for every "
                    "|omega|<1/2."
                ),
                blocker=None
                if original_closed
                else "Need original coordinate Weyl positivity.",
            ),
        },
        "originalWeylOperatorFamilyClosed": original_closed,
        "originalWeylOperatorPositiveClosed": original_closed,
        "originalWeylKernelPositivityClosed": original_closed,
        "proof": [
            "Import original coordinate Weyl kernel positivity for every |omega|<1/2.",
            "Import the kernel-to-operator identity in the fixed hbar=1 Weyl convention.",
            "Expose only this Weyl positive-operator family to the KLM packaging theorem.",
        ],
        "notClaimedHere": (
            "This theorem does not prove KLM positive type.  It only supplies "
            "the positive Weyl operator family consumed by the hbar=1 "
            "KLM/Weyl equivalence layer."
        ),
        "remainingAnalyticGap": None
        if original_closed
        else "The imported original Weyl positivity theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl positive-operator family theorem")
    print(f"  original Weyl input: {original_closed}")
    print(f"  operator family: {original_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
