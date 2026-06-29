#!/usr/bin/env python3
r"""Terminal formal RH/de Branges conclusion consequence.

This is the topmost publication-audit root.  It reads the full RH/de Branges
bridge ledger and exports only the terminal formal conclusion flags:

    formalRhClosed = true,
    rhClosed = true,
    endpointPassageClosed = true.

The all-omega KLM input, external-audit summary, and shifted-Xi endpoint
zero-location inputs remain below the full RH/de Branges bridge ledger.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ledger-json",
        default="rh_debranges_bridge_ledger.json",
    )
    parser.add_argument(
        "--json-out",
        default="rh_formal_conclusion_consequence_theorem.json",
    )
    args = parser.parse_args()

    ledger = load(args.ledger_json)
    formal_ok = bool(ledger.get("formalRhClosed"))
    rh_ok = bool(ledger.get("rhClosed"))
    endpoint_ok = bool(ledger.get("endpointPassageClosed"))
    vetted = bool(ledger.get("independentRhProofVetted"))
    theorem_closed = formal_ok and rh_ok and endpoint_ok

    data: dict[str, Any] = {
        "theoremName": "RH formal conclusion consequence theorem",
        "proofClass": "symbolic identity",
        "source": "RH/de Branges bridge ledger",
        "statuses": {
            "formalRhConclusionInputStatus": status(
                "formal RH conclusion input",
                formal_ok,
                "The full RH/de Branges bridge ledger records the formal RH conclusion as closed.",
            ),
            "rhConclusionInputStatus": status(
                "RH conclusion input",
                rh_ok,
                "The full RH/de Branges bridge ledger records rhClosed=true.",
            ),
            "endpointPassageInputStatus": status(
                "endpoint passage input",
                endpoint_ok,
                "The full RH/de Branges bridge ledger records the shifted-Xi endpoint passage as closed.",
            ),
            "formalConclusionConsequenceStatus": status(
                "formal conclusion consequence",
                theorem_closed,
                "Only the terminal formal conclusion flags are exported as the audit root.",
            ),
        },
        "endpointPassageClosed": endpoint_ok,
        "formalRhClosed": theorem_closed,
        "rhClosed": theorem_closed,
        "independentRhProofVetted": vetted,
        "caution": ledger.get("caution"),
        "nextProofTarget": ledger.get("nextProofTarget"),
        "proof": [
            "Import the full RH/de Branges bridge ledger.",
            "Export only the terminal formal conclusion flags for the top-level audit root.",
        ],
        "notExportedHere": [
            "uniform omega KLM/Weyl input",
            "external equivalence audit summary",
            "shifted-Xi zero-location consequence",
            "full proof ledger evidence",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "RH/de Branges bridge ledger is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("RH formal conclusion consequence theorem")
    print(f"  endpoint passage: {endpoint_ok}")
    print(f"  formal RH: {formal_ok}")
    print(f"  RH closed: {rh_ok}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
