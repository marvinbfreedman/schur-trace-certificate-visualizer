#!/usr/bin/env python3
r"""Narrow consequence of the augmented repair transport layer.

This file exposes only the consequence needed by
``xi_augmented_repair_positive_consequence_theorem.json``:

    D_aug >= 0 on X_aug
    K + R_aug^* D_aug R_aug is a closed nonnegative repaired form.

The detailed Mosco/closed-form transport proof remains in
``xi_augmented_repair_transport_limit_theorem.json`` and the completed tail
positivity proof remains in ``volterra_tail_positive_form_theorem.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


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
        "--transport-limit-json",
        default="xi_augmented_repair_transport_limit_theorem.json",
    )
    parser.add_argument(
        "--tail-positive-json",
        default="volterra_tail_positive_form_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_repair_transport_consequence_theorem.json",
    )
    args = parser.parse_args()

    limit = load(args.transport_limit_json)
    tail = load(args.tail_positive_json)

    limit_nonnegative_ok = bool(limit.get("DaugNonnegativeTransportLimitClosed"))
    limit_closed_form_ok = bool(limit.get("repairedAugmentedFormLimitClosed"))
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )
    transport_consequence_ok = limit_nonnegative_ok and limit_closed_form_ok and tail_ok

    data = {
        "theoremName": "augmented repair transport consequence theorem",
        "proofClass": "symbolic identity",
        "transportLimitJson": args.transport_limit_json,
        "tailPositiveJson": args.tail_positive_json,
        "transportedTraceSpace": {
            "name": "X_aug",
            "definition": "closure Ran(R_aug) with quotient norm inf{||f||_V:R_aug f=x}",
            "formDomain": "completed Volterra/Mellin graph-form domain V",
        },
        "statuses": {
            "transportLimitNonnegativeInputStatus": status(
                "transport-limit D_aug nonnegativity",
                limit_nonnegative_ok,
                "The transport-limit theorem gives a nonnegative trace-side limit form D_aug on X_aug.",
            ),
            "transportLimitClosedFormInputStatus": status(
                "closed repaired form limit",
                limit_closed_form_ok,
                "The transport-limit theorem gives the closed repaired augmented form as the finite repaired-form limit.",
            ),
            "positiveTailInputStatus": status(
                "positive Volterra tail input",
                tail_ok,
                "The positive Volterra tail theorem supplies the remaining completed-domain nonnegative tail input.",
            ),
            "transportConsequenceStatus": status(
                "transport consequence",
                transport_consequence_ok,
                "The two transport-limit conclusions and the positive tail input imply the narrow transport consequence.",
            ),
        },
        "transportedRepairClosed": transport_consequence_ok,
        "boundedNonnegativeRepairClosed": transport_consequence_ok,
        "closedRepairedAugmentedFormClosed": transport_consequence_ok,
        "DaugNonnegativeTransportLimitClosed": limit_nonnegative_ok,
        "DaugTransportConsequenceClosed": transport_consequence_ok,
        "proof": [
            "Import the closed-form transport-limit theorem.",
            "Import the positive Volterra tail theorem.",
            "Expose only D_aug>=0 on X_aug and the closed repaired-form conclusion.",
        ],
        "notExportedHere": [
            "finite Schur algebra",
            "Mosco graph-form proof",
            "trace quotient compatibility proof",
            "Volterra Green-lift proof",
        ],
        "remainingAnalyticGap": None
        if transport_consequence_ok
        else "Transport limit or positive tail input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented repair transport consequence theorem")
    print(f"  D_aug transport limit: {limit_nonnegative_ok}")
    print(f"  closed repaired form limit: {limit_closed_form_ok}")
    print(f"  positive tail: {tail_ok}")
    print(f"  consequence: {transport_consequence_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
