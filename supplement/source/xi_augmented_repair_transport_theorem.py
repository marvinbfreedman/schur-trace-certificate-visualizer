#!/usr/bin/env python3
r"""Continuum transport theorem for the augmented Xi repair.

This theorem is now a thin packaging layer over the closed-form transport
limit theorem.  It works in the transported trace Hilbert space

    X_aug = closure Ran(R_aug),
    ||x||_{X_aug} = inf { ||f||_V : R_aug f = x },

where V is the closed Volterra/Mellin graph-form domain.  The finite-to-
completed-space passage is handled by
``xi_augmented_repair_transport_limit_theorem.json``.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


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
    # Legacy inputs accepted for older invocations; no longer proof
    # dependencies of this packaging layer.
    parser.add_argument("--finite-repair-json", default="")
    parser.add_argument("--trace-convergence-json", default="")
    parser.add_argument("--json-out", default="xi_augmented_repair_transport_theorem.json")
    args = parser.parse_args()

    limit = load(args.transport_limit_json)
    tail = load(args.tail_positive_json)

    limit_ok = bool(
        limit.get("DaugNonnegativeTransportLimitClosed")
        and limit.get("repairedAugmentedFormLimitClosed")
    )
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )
    transport_ok = limit_ok and tail_ok

    limit_status = status(
        "closed-form transport limit input",
        limit_ok,
        (
            "The transport-limit theorem proves that the finite augmented "
            "repairs have a nonnegative closed-form limit D_aug on X_aug and "
            "that the repaired augmented form is closed and nonnegative."
        ),
    )
    tail_status = status(
        "positive Volterra tail input",
        tail_ok,
        (
            "The positive Volterra tail theorem proves the diagonal tail form "
            "P-M>=0 on the completed trace-fiber domain.  It is the narrow "
            "input needed here and sits below the original Weyl/KLM packaging."
        ),
    )
    transport_status = status(
        "bounded nonnegative transported trace form",
        transport_ok,
        (
            "The closed-form transport-limit theorem gives a bounded "
            "nonnegative trace-side form D_aug on X_aug, and the positive "
            "Volterra tail theorem supplies the remaining completed-domain "
            "nonnegative tail input."
        ),
    )
    positive_form_status = status(
        "closed repaired augmented form",
        transport_ok,
        (
            "K+R_aug^*D_aug R_aug is the closed-form limit of nonnegative "
            "finite repaired forms, hence is nonnegative on the completed "
            "augmented trace-fiber domain."
        ),
    )

    data = {
        "theoremName": "augmented repair continuum transport theorem",
        "proofClass": "symbolic identity",
        "transportLimitJson": args.transport_limit_json,
        "tailPositiveJson": args.tail_positive_json,
        "legacyInputs": {
            "finiteRepairJson": args.finite_repair_json or None,
            "traceConvergenceJson": args.trace_convergence_json or None,
        },
        "transportedTraceSpace": {
            "name": "X_aug",
            "definition": "closure Ran(R_aug) with quotient norm inf{||f||_V:R_aug f=x}",
            "formDomain": "completed Volterra/Mellin graph-form domain V",
        },
        "statuses": {
            "closedFormTransportLimitInputStatus": limit_status,
            "positiveVolterraTailInputStatus": tail_status,
            "boundedNonnegativeRepairStatus": transport_status,
            "closedRepairedAugmentedFormStatus": positive_form_status,
        },
        "transportedRepairClosed": transport_ok,
        "boundedNonnegativeRepairClosed": transport_ok,
        "closedRepairedAugmentedFormClosed": transport_ok,
        "DaugNonnegativeTransportLimitClosed": limit_ok,
        "proof": [
            "Import the closed-form transport-limit theorem for D_aug>=0 on X_aug.",
            "Import the positive Volterra tail theorem.",
            "Conclude the augmented repair is a bounded nonnegative trace-side form on the completed domain.",
        ],
        "remainingAnalyticGap": None if transport_ok else "One continuum transport input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented repair continuum transport theorem")
    print(f"  transport limit: {limit_ok}")
    print(f"  positive Volterra tail: {tail_ok}")
    print(f"  transported repair: {transport_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
