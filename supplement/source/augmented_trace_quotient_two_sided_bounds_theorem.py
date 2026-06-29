#!/usr/bin/env python3
r"""Two-sided augmented trace quotient bounds theorem."""

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
        default="augmented_trace_quotient_limsup_consequence_theorem.json",
    )
    parser.add_argument(
        "--quotient-liminf-json",
        default="augmented_trace_quotient_liminf_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_quotient_two_sided_bounds_theorem.json",
    )
    args = parser.parse_args()

    limsup = load(args.quotient_limsup_json)
    liminf = load(args.quotient_liminf_json)

    upper_ok = bool(limsup.get("traceQuotientLimsupClosed"))
    lower_ok = bool(liminf.get("traceQuotientLiminfClosed"))
    two_sided_ok = upper_ok and lower_ok

    data = {
        "theoremName": "augmented trace quotient two-sided bounds theorem",
        "proofClass": "analytic proof",
        "quotientLimsupJson": args.quotient_limsup_json,
        "quotientLiminfJson": args.quotient_liminf_json,
        "statement": (
            "The augmented trace quotient norms satisfy both the upper recovery "
            "bound and the lower representative-compactness bound."
        ),
        "statuses": {
            "quotientUpperBoundInputStatus": status(
                "trace quotient upper bound input",
                upper_ok,
                "The quotient limsup theorem supplies the recovery upper bound.",
            ),
            "quotientLowerBoundInputStatus": status(
                "trace quotient lower bound input",
                lower_ok,
                "The quotient liminf theorem supplies the representative lower bound.",
            ),
            "twoSidedQuotientBoundsStatus": status(
                "two-sided trace quotient bounds",
                two_sided_ok,
                (
                    "The upper and lower quotient bounds identify the limiting "
                    "transported trace quotient norm from both sides."
                ),
            ),
        },
        "twoSidedTraceQuotientBoundsClosed": two_sided_ok,
        "traceQuotientLimsupClosed": upper_ok,
        "traceQuotientLiminfClosed": lower_ok,
        "proof": [
            "Import the trace quotient limsup theorem.",
            "Import the trace quotient liminf theorem.",
            "Use the two inequalities as the two-sided quotient bound package.",
        ],
        "remainingAnalyticGap": None if two_sided_ok else "A quotient bound input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient two-sided bounds theorem")
    print(f"  quotient upper bound: {upper_ok}")
    print(f"  quotient lower bound: {lower_ok}")
    print(f"  two-sided bounds: {two_sided_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
