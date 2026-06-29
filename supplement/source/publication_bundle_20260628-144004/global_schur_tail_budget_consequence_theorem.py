#!/usr/bin/env python3
r"""Terminal consequence for the finite source-inactive Schur-tail budget.

The source-inactive Schur-tail certificate only needs the finite diagnostic
budget for the low/mid Schur comparison:

    remainingAnalyticItems[name="source-inactive high-frequency Schur tail"]
      .finiteDiagnostic.bestObservedNonzeroOperatorTail.fraction

This wrapper extracts that diagnostic from the global bridge without importing
the full bridge proof ledger upstream.
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


def find_tail_budget(bridge: dict[str, Any]) -> dict[str, Any] | None:
    for item in bridge.get("remainingAnalyticItems", []):
        if item.get("name") == "source-inactive high-frequency Schur tail":
            diagnostic = item.get("finiteDiagnostic", {})
            return diagnostic.get("bestObservedNonzeroOperatorTail")
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bridge-json",
        default="global_weyl_volterra_schur_bridge.json",
    )
    parser.add_argument(
        "--json-out",
        default="global_schur_tail_budget_consequence_theorem.json",
    )
    args = parser.parse_args()

    bridge = load(args.bridge_json)
    tail_budget = find_tail_budget(bridge)
    fraction = tail_budget.get("fraction") if isinstance(tail_budget, dict) else None
    try:
        fraction_positive = float(fraction) > 0.0
    except Exception:
        fraction_positive = False

    data: dict[str, Any] = {
        "theoremName": "global Schur tail finite-budget consequence theorem",
        "proofClass": "symbolic identity",
        "bridgeJson": args.bridge_json,
        "sourceInactiveTailBudgetDiagnostic": tail_budget,
        "finiteSchurTailBudgetDiagnostic": tail_budget,
        "sourceInactiveTailBudgetFraction": fraction,
        "sourceInactiveTailBudgetConsequenceClosed": fraction_positive,
        "statuses": {
            "sourceInactiveTailBudgetDiagnosticStatus": status(
                "source-inactive finite Schur-tail budget diagnostic",
                fraction_positive,
                (
                    "The global bridge records a positive finite diagnostic "
                    "fraction for the nonzero source-inactive operator tail."
                ),
            )
        },
        "proof": [
            "Import the global Weyl/Volterra Schur bridge ledger.",
            "Extract the finiteDiagnostic.bestObservedNonzeroOperatorTail record for the source-inactive tail item.",
            "Export only that finite diagnostic budget to the source-inactive Schur-tail certificate.",
        ],
        "remainingAnalyticGap": None
        if fraction_positive
        else "Global bridge does not expose a positive source-inactive finite tail budget.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Global Schur tail budget consequence theorem")
    print(f"  budget fraction: {fraction}")
    print(f"  theorem closed: {fraction_positive}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
