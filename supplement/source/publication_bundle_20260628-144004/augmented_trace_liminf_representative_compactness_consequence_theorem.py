#!/usr/bin/env python3
r"""Narrow consequence of augmented trace liminf representative compactness."""

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
        "--representative-compactness-json",
        default="augmented_trace_liminf_representative_compactness_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_liminf_representative_compactness_consequence_theorem.json",
    )
    args = parser.parse_args()

    representative = load(args.representative_compactness_json)
    representative_ok = bool(
        representative.get("representativeCompactnessAndTraceIdentificationClosed")
    )
    weak_ok = bool(representative.get("weakRepresentativeCompactnessClosed", representative_ok))
    trace_ok = bool(representative.get("traceLimitIdentificationClosed", representative_ok))
    ok = representative_ok and weak_ok and trace_ok

    data = {
        "theoremName": "augmented trace liminf representative compactness consequence theorem",
        "proofClass": "symbolic identity",
        "representativeCompactnessJson": args.representative_compactness_json,
        "statement": (
            "The augmented trace liminf representative compactness theorem "
            "supplies the representative compactness and trace-identification "
            "consequence needed by the quotient liminf theorem."
        ),
        "statuses": {
            "weakRepresentativeCompactnessStatus": status(
                "weak representative compactness",
                weak_ok,
                "Imported from the augmented trace liminf representative compactness theorem.",
            ),
            "traceLimitIdentificationStatus": status(
                "trace limit identification",
                trace_ok,
                "Imported from the augmented trace liminf representative compactness theorem.",
            ),
            "representativeCompactnessAndTraceIdentificationStatus": status(
                "representative compactness and trace identification",
                representative_ok,
                "Imported from the augmented trace liminf representative compactness theorem.",
            ),
            "representativeCompactnessConsequenceStatus": status(
                "representative compactness consequence",
                ok,
                (
                    "Bounded near-minimizing finite representatives have weak "
                    "graph limits whose augmented traces identify the limiting "
                    "trace coordinate."
                ),
            ),
        },
        "representativeCompactnessAndTraceIdentificationClosed": representative_ok,
        "weakRepresentativeCompactnessClosed": weak_ok,
        "traceLimitIdentificationClosed": trace_ok,
        "representativeCompactnessConsequenceClosed": ok,
        "proof": [
            "Import the augmented trace liminf representative compactness theorem.",
            "Extract the representative compactness and trace-identification conclusion.",
            "Expose only that consequence to the augmented trace quotient liminf theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Representative compactness or trace-limit identification is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace liminf representative compactness consequence theorem")
    print(f"  weak representative compactness: {weak_ok}")
    print(f"  trace limit identification: {trace_ok}")
    print(f"  representative consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
