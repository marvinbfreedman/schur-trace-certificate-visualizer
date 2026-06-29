#!/usr/bin/env python3
r"""Terminal consequence of the global Weyl/Volterra Schur bridge.

The external equivalence audit only needs the closed conclusion that the
global Weyl/Volterra Schur theorem has been assembled from the full-Phi source
noncollapse, active-range, source-inactive tail, high-block exhaustion, and
quotient Schur inputs.  The detailed bridge ledger remains one layer lower.
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
    parser.add_argument(
        "--bridge-json",
        default="global_weyl_volterra_schur_bridge.json",
    )
    parser.add_argument(
        "--json-out",
        default="global_weyl_volterra_schur_bridge_consequence_theorem.json",
    )
    args = parser.parse_args()

    bridge = load(args.bridge_json)
    closed = bool(bridge.get("globalSchurTheoremClosed"))

    bridge_status = status(
        "global Weyl/Volterra Schur bridge consequence",
        closed,
        (
            "The global bridge ledger closes the full-Phi source noncollapse, "
            "active range inclusion, source-inactive tail absorption, "
            "high-block exhaustion, and quotient Schur assembly needed for "
            "the normalized Weyl/Volterra Schur certificate."
        ),
        blocker=None if closed else "Close global_weyl_volterra_schur_bridge.py.",
    )

    data = {
        "theoremName": "global Weyl/Volterra Schur bridge consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "globalSchurBridgeConsequenceStatus": bridge_status,
        },
        "globalSchurTheoremClosed": bridge_status["closed"],
        "globalWeylVolterraSchurBridgeClosed": bridge_status["closed"],
        "bridgePasses": bool(bridge.get("bridgePasses")),
        "proof": [
            "Import the detailed global Weyl/Volterra Schur bridge ledger.",
            "Export only the terminal Schur bridge consequence needed by the external equivalence audit.",
        ],
        "remainingAnalyticGap": None if bridge_status["closed"] else bridge_status["blocker"],
        "nextProofTarget": None if bridge_status["closed"] else bridge_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Global Weyl/Volterra Schur bridge consequence theorem")
    print(f"  global Schur theorem closed: {data['globalSchurTheoremClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
