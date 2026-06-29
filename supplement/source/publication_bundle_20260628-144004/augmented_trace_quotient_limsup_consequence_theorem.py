#!/usr/bin/env python3
r"""Narrow consequence of the augmented trace quotient limsup theorem."""

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
        "--quotient-limsup-json",
        default="augmented_trace_quotient_limsup_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_quotient_limsup_consequence_theorem.json",
    )
    args = parser.parse_args()

    limsup = load(args.quotient_limsup_json)
    limsup_ok = bool(limsup.get("traceQuotientLimsupClosed"))
    upper_ok = bool(limsup.get("finiteQuotientUpperBoundClosed", limsup_ok))
    ok = limsup_ok and upper_ok

    data = {
        "theoremName": "augmented trace quotient limsup consequence theorem",
        "proofClass": "symbolic identity",
        "quotientLimsupJson": args.quotient_limsup_json,
        "statement": (
            "The augmented trace quotient limsup theorem supplies the quotient "
            "upper-bound consequence needed by the two-sided quotient bounds theorem."
        ),
        "statuses": {
            "traceQuotientLimsupStatus": status(
                "trace quotient limsup",
                limsup_ok,
                "Imported from the augmented trace quotient limsup theorem.",
            ),
            "finiteQuotientUpperBoundStatus": status(
                "finite quotient upper bound",
                upper_ok,
                "Imported from the augmented trace quotient limsup theorem.",
            ),
            "traceQuotientLimsupConsequenceStatus": status(
                "trace quotient limsup consequence",
                ok,
                "The transported trace quotient upper bound is available upstream.",
            ),
        },
        "traceQuotientLimsupClosed": limsup_ok,
        "finiteQuotientUpperBoundClosed": upper_ok,
        "traceQuotientLimsupConsequenceClosed": ok,
        "proof": [
            "Import the augmented trace quotient limsup theorem.",
            "Extract the trace quotient upper-bound conclusion.",
            "Expose only that consequence to the two-sided quotient bounds theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Trace quotient limsup input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient limsup consequence theorem")
    print(f"  trace quotient limsup: {limsup_ok}")
    print(f"  quotient upper bound: {upper_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
