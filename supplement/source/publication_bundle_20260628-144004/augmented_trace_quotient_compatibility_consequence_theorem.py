#!/usr/bin/env python3
r"""Narrow consequence of augmented trace quotient compatibility."""

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
        "--trace-quotient-json",
        default="augmented_trace_quotient_compatibility_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_quotient_compatibility_consequence_theorem.json",
    )
    args = parser.parse_args()

    quotient = load(args.trace_quotient_json)
    quotient_ok = bool(quotient.get("traceQuotientCompatibilityClosed"))
    norm_ok = bool(quotient.get("finiteTraceQuotientNormConvergenceClosed", quotient_ok))
    ok = quotient_ok and norm_ok

    data = {
        "theoremName": "augmented trace quotient compatibility consequence theorem",
        "proofClass": "symbolic identity",
        "traceQuotientJson": args.trace_quotient_json,
        "statement": (
            "The augmented trace quotient compatibility theorem supplies the "
            "transported quotient compatibility consequence needed by null-fiber "
            "repair compatibility."
        ),
        "statuses": {
            "traceQuotientCompatibilityStatus": status(
                "trace quotient compatibility",
                quotient_ok,
                "Imported from the augmented trace quotient compatibility theorem.",
            ),
            "finiteTraceQuotientNormConvergenceStatus": status(
                "finite trace quotient norm convergence",
                norm_ok,
                "Imported from the augmented trace quotient compatibility theorem.",
            ),
            "traceQuotientCompatibilityConsequenceStatus": status(
                "trace quotient compatibility consequence",
                ok,
                (
                    "Finite trace equivalence is identified with the completed "
                    "transported quotient space X_aug."
                ),
            ),
        },
        "traceQuotientCompatibilityClosed": quotient_ok,
        "finiteTraceQuotientNormConvergenceClosed": norm_ok,
        "traceQuotientCompatibilityConsequenceClosed": ok,
        "proof": [
            "Import the augmented trace quotient compatibility theorem.",
            "Extract the transported trace quotient compatibility conclusion.",
            "Expose only that consequence to null-fiber repair compatibility.",
        ],
        "remainingAnalyticGap": None if ok else "Trace quotient compatibility input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient compatibility consequence theorem")
    print(f"  trace quotient compatibility: {quotient_ok}")
    print(f"  quotient norm convergence: {norm_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
