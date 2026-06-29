#!/usr/bin/env python3
r"""Terminal consequence of the RH/de Branges bridge ledger.

The external Weyl/Volterra audit only needs the final formal bridge flag from
the RH/de Branges ledger.  The uniform-omega KLM input, external audit
consequence, shifted-Xi kernel positivity, and endpoint zero-location inputs
remain audited one layer lower.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rh-ledger-json", default="rh_debranges_bridge_ledger.json")
    parser.add_argument(
        "--json-out",
        default="rh_debranges_bridge_ledger_consequence_theorem.json",
    )
    args = parser.parse_args()

    ledger = load(args.rh_ledger_json)
    formal_closed = bool(ledger.get("formalRhClosed"))

    rh_status = status(
        "RH/de Branges bridge ledger consequence",
        formal_closed,
        (
            "The RH/de Branges bridge ledger records the formal implication "
            "from the all-omega KLM/Weyl positive-type input and shifted-Xi "
            "endpoint zero-location consequence to the RH-side conclusion."
        ),
        blocker=None if formal_closed else "Close rh_debranges_bridge_ledger.py.",
    )

    data = {
        "theoremName": "RH/de Branges bridge ledger consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "rhBridgeLedgerConsequenceStatus": rh_status,
        },
        "formalRhClosed": rh_status["closed"],
        "conditionalRhClosed": bool(ledger.get("conditionalRhClosed") or formal_closed),
        "rhBridgeLedgerClosed": rh_status["closed"],
        "proof": [
            "Import the detailed RH/de Branges bridge ledger.",
            "Export only the terminal formal RH bridge consequence needed by the external audit.",
        ],
        "remainingAnalyticGap": None if rh_status["closed"] else rh_status["blocker"],
        "nextProofTarget": None if rh_status["closed"] else rh_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("RH/de Branges bridge ledger consequence theorem")
    print(f"  formal RH bridge closed: {data['formalRhClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
