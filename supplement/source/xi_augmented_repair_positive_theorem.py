#!/usr/bin/env python3
r"""Legacy wrapper for augmented Xi-trace repair positivity.

The narrow theorem object used upstream is now
``xi_augmented_repair_positive_consequence_theorem.json``.  This compatibility
wrapper keeps the older filename available while avoiding a direct proof edge
from the closed-cone bridge to the detailed transport machinery.
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
        "--consequence-json",
        default="xi_augmented_repair_positive_consequence_theorem.json",
    )
    # Legacy argument accepted for old invocations.  It is no longer a proof
    # dependency unless explicitly supplied together with a regenerated
    # consequence theorem.
    parser.add_argument("--transport-json", default="")
    parser.add_argument("--json-out", default="xi_augmented_repair_positive_theorem.json")
    args = parser.parse_args()

    consequence = load(args.consequence_json)

    repair_ok = bool(consequence.get("continuumAugmentedRepairClosed"))
    nonnegative_ok = bool(
        consequence.get("DaugNonnegativeOnTransportedTraceSpaceClosed")
        or consequence.get("DaugBoundedNonnegativeClosed")
    )
    closed_form_ok = bool(consequence.get("closedAugmentedPositiveFormClosed"))

    nonnegative_status = status(
        "D_aug nonnegative consequence",
        nonnegative_ok,
        (
            "Imported from the narrow consequence theorem: D_aug is bounded "
            "and nonnegative on the completed transported augmented trace "
            "space."
        ),
    )
    closed_form_status = status(
        "closed augmented positive form",
        closed_form_ok,
        (
            "Imported from the narrow consequence theorem: the repaired "
            "augmented form is closed on the completed trace-fiber domain."
        ),
    )
    wrapper_status = status(
        "legacy repair positivity wrapper",
        repair_ok,
        (
            "This file is a backward-compatible alias for the narrow positive "
            "repair consequence and no longer exposes the internal transport "
            "proof as a direct upstream dependency."
        ),
    )

    data = {
        "theoremName": "augmented Xi-trace repair positivity theorem legacy wrapper",
        "proofClass": "analytic proof wrapper",
        "repairPositiveConsequenceJson": args.consequence_json,
        "legacyInputs": {
            "transportJson": args.transport_json or None,
        },
        "transportedTraceSpace": consequence.get("transportedTraceSpace"),
        "statuses": {
            "DaugNonnegativeConsequenceStatus": nonnegative_status,
            "closedAugmentedPositiveFormStatus": closed_form_status,
            "legacyRepairPositivityWrapperStatus": wrapper_status,
        },
        "repairPositivityClosed": repair_ok,
        "continuumAugmentedRepairClosed": repair_ok,
        "DaugNonnegativeOnTransportedTraceSpaceClosed": nonnegative_ok,
        "closedAugmentedPositiveFormClosed": closed_form_ok,
        "proof": [
            "Import the narrow positive-repair consequence theorem.",
            "Re-export legacy compatibility flags expected by older bridge scripts.",
        ],
        "remainingAnalyticGap": None
        if repair_ok
        else "The narrow positive-repair consequence theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented Xi-trace repair positivity theorem")
    print(f"  D_aug nonnegative consequence: {nonnegative_ok}")
    print(f"  closed repaired form: {closed_form_ok}")
    print(f"  repair positivity: {repair_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
