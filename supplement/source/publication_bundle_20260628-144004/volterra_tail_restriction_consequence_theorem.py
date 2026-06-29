#!/usr/bin/env python3
r"""Fixed-core restriction consequence of Volterra tail positivity.

This wrapper exposes only the input needed by the augmented pointwise repaired
form certificate:

    t_N = T_volt restricted to the fixed augmented trace core V_N is >= 0.

The full Volterra tail positivity proof remains in
``volterra_tail_positive_form_theorem.json``.
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
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tail-positive-json",
        default="volterra_tail_positive_form_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="volterra_tail_restriction_consequence_theorem.json",
    )
    args = parser.parse_args()

    tail = load(args.tail_positive_json)
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )
    ok = tail_ok

    data = {
        "theoremName": "Volterra tail fixed-core restriction consequence",
        "proofClass": "symbolic identity",
        "tailPositiveJson": args.tail_positive_json,
        "statement": (
            "Restricting the nonnegative Volterra tail quadratic form to any "
            "fixed augmented trace core preserves nonnegativity."
        ),
        "statuses": {
            "volterraTailPositiveInputStatus": status(
                "Volterra tail positivity input",
                tail_ok,
                "The full Volterra tail theorem proves T_volt=P-M>=0.",
            ),
            "fixedCoreRestrictionStatus": status(
                "fixed-core tail restriction",
                ok,
                (
                    "If a quadratic form is nonnegative on a completed domain, "
                    "then its restriction to any fixed finite/core subspace is "
                    "nonnegative."
                ),
            ),
        },
        "volterraTailRestrictionNonnegativeClosed": ok,
        "fixedCoreOnly": True,
        "proof": [
            "Import the Volterra tail positive form theorem.",
            "Let V_N be any fixed augmented trace core contained in the form domain.",
            "For f in V_N, t_N(f)=T_volt(f)>=0 by the full-domain positivity.",
        ],
        "remainingAnalyticGap": None if ok else "Volterra tail positivity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Volterra tail restriction consequence theorem")
    print(f"  tail positivity input: {tail_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
