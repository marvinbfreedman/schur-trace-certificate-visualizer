#!/usr/bin/env python3
r"""Trace quotient compatibility theorem for augmented transport."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--two-sided-bounds-json",
        default="augmented_trace_quotient_two_sided_bounds_theorem.json",
    )
    parser.add_argument(
        "--quotient-limsup-json",
        default="",
    )
    parser.add_argument(
        "--quotient-liminf-json",
        default="",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_quotient_compatibility_theorem.json",
    )
    args = parser.parse_args()

    bounds = load(args.two_sided_bounds_json)

    limsup_ok = bool(bounds.get("traceQuotientLimsupClosed"))
    liminf_ok = bool(bounds.get("traceQuotientLiminfClosed"))
    bounds_ok = bool(bounds.get("twoSidedTraceQuotientBoundsClosed"))
    quotient_ok = bounds_ok and limsup_ok and liminf_ok

    data = {
        "theoremName": "augmented trace quotient compatibility identity",
        "proofClass": "symbolic identity",
        "twoSidedBoundsJson": args.two_sided_bounds_json,
        "legacyQuotientLimsupJson": args.quotient_limsup_json or None,
        "legacyQuotientLiminfJson": args.quotient_liminf_json or None,
        "transportedTraceSpace": {
            "name": "X_aug",
            "definition": "closure Ran(R_aug) with quotient norm inf{||f||_V:R_aug f=x}",
        },
        "statement": (
            "The finite trace quotient norms induced by R_aug,N converge to "
            "the transported quotient norm on X_aug in the Mosco sense."
        ),
        "statuses": {
            "twoSidedBoundsInputStatus": status(
                "two-sided trace quotient bounds input",
                bounds_ok,
                (
                    "The two-sided bounds theorem packages the upper and "
                    "lower quotient inequalities."
                ),
            ),
            "quotientLimsupInputStatus": status(
                "trace quotient upper bound input",
                limsup_ok,
                "The two-sided bounds theorem includes the recovery upper bound.",
            ),
            "quotientLiminfInputStatus": status(
                "trace quotient lower bound input",
                liminf_ok,
                "The two-sided bounds theorem includes the representative compactness lower bound.",
            ),
            "traceQuotientCompatibilityStatus": status(
                "trace quotient compatibility",
                quotient_ok,
                (
                    "By definition, having both the upper and lower quotient "
                    "bounds is exactly trace quotient compatibility on X_aug."
                ),
            ),
        },
        "traceQuotientCompatibilityClosed": quotient_ok,
        "finiteTraceQuotientNormConvergenceClosed": quotient_ok,
        "proof": [
            "Import the two-sided quotient bounds theorem.",
            "Read off the upper and lower quotient inequalities.",
            "Use the definition of trace quotient compatibility.",
        ],
        "remainingAnalyticGap": None if quotient_ok else "A quotient compatibility input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient compatibility theorem")
    print(f"  quotient limsup: {limsup_ok}")
    print(f"  quotient liminf: {liminf_ok}")
    print(f"  quotient compatibility: {quotient_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
