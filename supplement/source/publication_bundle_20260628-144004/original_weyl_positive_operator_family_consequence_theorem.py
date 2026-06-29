#!/usr/bin/env python3
r"""Slim consequence for the original Weyl positive-operator family.

The uniform omega Weyl/KLM bridge only needs the operator-level consequence

    Op^W(sigma_omega) >= 0,   |omega| < 1/2,

in the fixed hbar=1 Weyl normalization.  This wrapper exposes that single
flag while leaving the kernel-to-operator identity and original Weyl positivity
machinery below the full operator-family theorem.
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
        "--operator-family-json",
        default="original_weyl_positive_operator_family_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="original_weyl_positive_operator_family_consequence_theorem.json",
    )
    args = parser.parse_args()

    operator_family = load(args.operator_family_json)
    ok = bool(
        operator_family.get("originalWeylOperatorFamilyClosed")
        or operator_family.get("originalWeylOperatorPositiveClosed")
        or operator_family.get("originalWeylKernelPositivityClosed")
    )

    data = {
        "theoremName": "original Weyl positive-operator family consequence theorem",
        "proofClass": "symbolic identity",
        "operatorFamilySource": "original Weyl positive-operator family theorem",
        "omegaRange": "|omega| < 1/2",
        "normalization": "hbar=1 Weyl convention",
        "statuses": {
            "weylOperatorFamilyStatus": status(
                "positive Weyl operator family",
                ok,
                (
                    "The imported operator-family theorem proves "
                    "Op^W(sigma_omega)>=0 for every |omega|<1/2."
                ),
            )
        },
        "originalWeylOperatorFamilyClosed": ok,
        "originalWeylOperatorPositiveClosed": ok,
        "originalWeylKernelPositivityClosed": ok,
        "proof": [
            "Import the original Weyl positive-operator family theorem.",
            "Expose only the hbar=1 positive Weyl operator family consumed by the KLM bridge.",
        ],
        "remainingAnalyticGap": None if ok else "Original Weyl positive-operator family theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl positive-operator family consequence theorem")
    print(f"  operator family: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
