#!/usr/bin/env python3
r"""Narrow consequence of the augmented trace quotient liminf theorem."""

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
        "--quotient-liminf-json",
        default="augmented_trace_quotient_liminf_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_quotient_liminf_consequence_theorem.json",
    )
    args = parser.parse_args()

    liminf = load(args.quotient_liminf_json)
    liminf_ok = bool(liminf.get("traceQuotientLiminfClosed"))
    lower_ok = bool(liminf.get("finiteQuotientLowerBoundClosed", liminf_ok))
    ok = liminf_ok and lower_ok

    data = {
        "theoremName": "augmented trace quotient liminf consequence theorem",
        "proofClass": "symbolic identity",
        "quotientLiminfJson": args.quotient_liminf_json,
        "statement": (
            "The augmented trace quotient liminf theorem supplies the quotient "
            "lower-bound consequence needed by the two-sided quotient bounds theorem."
        ),
        "statuses": {
            "traceQuotientLiminfStatus": status(
                "trace quotient liminf",
                liminf_ok,
                "Imported from the augmented trace quotient liminf theorem.",
            ),
            "finiteQuotientLowerBoundStatus": status(
                "finite quotient lower bound",
                lower_ok,
                "Imported from the augmented trace quotient liminf theorem.",
            ),
            "traceQuotientLiminfConsequenceStatus": status(
                "trace quotient liminf consequence",
                ok,
                "The transported trace quotient lower bound is available upstream.",
            ),
        },
        "traceQuotientLiminfClosed": liminf_ok,
        "finiteQuotientLowerBoundClosed": lower_ok,
        "traceQuotientLiminfConsequenceClosed": ok,
        "proof": [
            "Import the augmented trace quotient liminf theorem.",
            "Extract the trace quotient lower-bound conclusion.",
            "Expose only that consequence to the two-sided quotient bounds theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Trace quotient liminf input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient liminf consequence theorem")
    print(f"  trace quotient liminf: {liminf_ok}")
    print(f"  quotient lower bound: {lower_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
