#!/usr/bin/env python3
r"""Uniform omega bridge from Weyl operator positivity to KLM positivity.

The original Weyl proof is consumed through the narrow theorem object
``original_weyl_positive_operator_family_theorem.json``.  That file exposes
only the positive Weyl operator family

    Op^W(sigma_omega) >= 0,  |omega|<1/2,

and keeps the Volterra/Green, primitive endpoint, and quotient-to-original lift
machinery one layer lower.  This bridge only applies the hbar=1 KLM/Weyl
equivalence theorem and packages the result as KLM quantum positive type for
Q_omega.

It deliberately does not claim the final RH/de Branges bridge; that remains a
separate theorem.
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
        "--weyl-positive-json",
        default="original_weyl_positive_operator_family_consequence_theorem.json",
    )
    # Legacy name accepted for old command lines; it is not used unless the
    # new narrow input is absent.
    parser.add_argument("--original-weyl-json", default="")
    parser.add_argument(
        "--klm-weyl-json",
        default="klm_weyl_hbar1_equivalence_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="uniform_omega_weyl_klm_bridge.json")
    args = parser.parse_args()

    weyl_positive_path = args.weyl_positive_json or args.original_weyl_json
    original = load(weyl_positive_path)
    klm_weyl = load(args.klm_weyl_json)

    weyl_operator_closed = bool(
        original.get("originalWeylOperatorPositiveClosed")
        or original.get("originalWeylOperatorFamilyClosed")
        or original.get("originalWeylKernelPositivityClosed")
        or nested_closed(original, "statuses", "weylOperatorFamilyStatus")
    )
    klm_normalization = bool(klm_weyl.get("klmWeylHbar1EquivalenceClosed"))

    original_status = status(
        "positive Weyl operator family for |omega|<1/2",
        weyl_operator_closed,
        (
            "Imported from the narrow original Weyl positive-operator family "
            "theorem: Op^W(sigma_omega) is positive semidefinite for every "
            "|omega|<1/2."
        ),
        blocker=None
        if weyl_operator_closed
        else "Close original_weyl_positive_operator_family_theorem.py.",
    )
    klm_status = status(
        "KLM quantum positive-type condition for Q_omega",
        weyl_operator_closed and klm_normalization,
        (
            "In the hbar=1 convention, the standalone KLM/Weyl equivalence "
            "theorem identifies original Weyl positivity with KLM quantum "
            "positive type of Q_omega."
        ),
        blocker=None
        if weyl_operator_closed and klm_normalization
        else "Need Weyl operator positivity and the hbar=1 KLM/Weyl equivalence.",
    )
    data = {
        "theoremName": "uniform omega Weyl/KLM bridge",
        "omegaRange": "|omega| < 1/2",
        "weylPositiveJson": weyl_positive_path,
        "legacyOriginalWeylJson": args.original_weyl_json or None,
        "klmWeylJson": args.klm_weyl_json,
        "statuses": {
            "originalWeylPositivityStatus": original_status,
            "klmPositiveTypeStatus": klm_status,
        },
        "proof": [
            "Import the positive Weyl operator family for every |omega|<1/2.",
            "Import the hbar=1 KLM/Weyl equivalence theorem.",
            "Package Op^W(sigma_omega)>=0 as KLM quantum positive type of Q_omega.",
        ],
        "uniformOmegaCoverageClosed": original_status["closed"],
        "originalWeylKernelPositivityClosed": original_status["closed"],
        "originalWeylOperatorPositiveClosed": original_status["closed"],
        "originalKlmConditionClosed": klm_status["closed"],
        "notClaimedHere": (
            "This theorem stops at KLM quantum positive type.  The RH-facing "
            "positive-kernel and endpoint consequences are handled by the top "
            "bridge ledger."
        ),
        "remainingAnalyticGap": None if klm_status["closed"] else klm_status["blocker"],
        "nextProofTarget": "Use the top bridge ledger for the downstream RH-facing layer.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Uniform omega Weyl/KLM bridge")
    print(f"  uniform omega coverage: {data['uniformOmegaCoverageClosed']}")
    print(f"  Weyl operator positivity: {original_status['closed']}")
    print(f"  KLM positive type: {klm_status['closed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
