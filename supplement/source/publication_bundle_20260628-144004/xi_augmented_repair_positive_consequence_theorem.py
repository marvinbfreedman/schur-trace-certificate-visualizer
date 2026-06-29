#!/usr/bin/env python3
r"""Narrow positive-repair consequence for the augmented Xi trace space.

This theorem is intentionally smaller than the transported-repair theorem it
imports.  It exposes only the consequence needed by the closed-cone bridge:

    D_aug >= 0 on X_aug = closure Ran(R_aug)

where X_aug carries the transported quotient trace norm.  The finite Schur
construction, trace convergence, Volterra tail positivity, and closed-form
transport machinery remain in ``xi_augmented_repair_transport_theorem.json``.
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
        "--transport-json",
        default="xi_augmented_repair_transport_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_repair_positive_consequence_theorem.json",
    )
    args = parser.parse_args()

    transport = load(args.transport_json)

    nonnegative_repair_ok = bool(
        (
            transport.get("boundedNonnegativeRepairClosed")
            or transport.get("DaugNonnegativeTransportLimitClosed")
        )
        and transport.get("transportedRepairClosed")
    )
    closed_form_ok = bool(transport.get("closedRepairedAugmentedFormClosed"))
    consequence_ok = nonnegative_repair_ok and closed_form_ok

    trace_space = transport.get("transportedTraceSpace") or {
        "name": "X_aug",
        "definition": "closure Ran(R_aug) with transported quotient norm",
    }

    nonnegative_status = status(
        "D_aug nonnegative on transported augmented trace space",
        nonnegative_repair_ok,
        (
            "The continuum transport theorem constructs a bounded nonnegative "
            "trace-side form D_aug on X_aug=closure Ran(R_aug), equipped with "
            "the quotient trace norm."
        ),
    )
    closed_form_status = status(
        "closed repaired augmented form",
        closed_form_ok,
        (
            "The repaired form K+R_aug^*D_aug R_aug is the closed lower-"
            "semicontinuous limit of nonnegative finite repaired forms."
        ),
    )
    consequence_status = status(
        "positive repair consequence used by closed cone",
        consequence_ok,
        (
            "Only the narrow consequence D_aug>=0 on X_aug, together with the "
            "closed repaired form, is exported to the de Branges closed-cone "
            "layer."
        ),
    )

    data = {
        "theoremName": "augmented Xi repair positive consequence theorem",
        "proofClass": "symbolic identity",
        "transportJson": args.transport_json,
        "transportedTraceSpace": trace_space,
        "exportedConsequence": {
            "operator": "D_aug",
            "space": "X_aug = closure Ran(R_aug)",
            "norm": "transported quotient trace norm",
            "claim": "D_aug >= 0 as a bounded nonnegative trace-side form",
        },
        "statuses": {
            "DaugNonnegativeOnTransportedTraceSpaceStatus": nonnegative_status,
            "closedRepairedAugmentedFormStatus": closed_form_status,
            "positiveRepairConsequenceStatus": consequence_status,
        },
        "DaugNonnegativeOnTransportedTraceSpaceClosed": nonnegative_repair_ok,
        "DaugBoundedNonnegativeClosed": nonnegative_repair_ok,
        "augmentedRepairPositiveConsequenceClosed": consequence_ok,
        "continuumAugmentedRepairClosed": consequence_ok,
        "closedAugmentedPositiveFormClosed": closed_form_ok,
        "proof": [
            "Import the continuum transported-repair theorem.",
            "Project away its construction details and retain only D_aug>=0 on X_aug.",
            "Retain the closed repaired-form conclusion needed for cone closure.",
        ],
        "notExportedHere": [
            "finite augmented Schur construction",
            "trace convergence proof",
            "Volterra tail positivity proof",
            "Galerkin transport machinery",
        ],
        "remainingAnalyticGap": None
        if consequence_ok
        else "The transported repair theorem has not closed the positive consequence.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented Xi repair positive consequence theorem")
    print(f"  D_aug nonnegative on X_aug: {nonnegative_repair_ok}")
    print(f"  closed repaired form: {closed_form_ok}")
    print(f"  exported consequence: {consequence_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
