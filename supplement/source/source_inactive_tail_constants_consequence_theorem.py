#!/usr/bin/env python3
r"""Terminal constants consequence for the source-inactive tail estimate.

The detailed interval certificate proves the finite source-model constants and
the low/mid Schur absorption comparison.  The high-block exhaustion theorem
only needs the resulting constants:

    normalizedEpsilonDelta < finiteLowMidSchurBudget,
    sourceInactiveTailConstantsClosed = true.

This wrapper exports only that consequence.  It is deliberately terminal in the
publication dependency graph: the raw interval certificates remain available as
audit artifacts, while the continuum high-block proof imports this exact input
contract.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {
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
        "--constants-json",
        default="source_inactive_tail_constants_theorem.json",
        help="Detailed constants theorem used to generate this terminal input.",
    )
    parser.add_argument(
        "--json-out",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    args = parser.parse_args()

    constants = load(args.constants_json)

    epsilon = float(constants["normalizedEpsilonDelta"])
    budget = float(constants["finiteLowMidSchurBudget"])
    slack = float(constants["absorptionSlack"])
    constants_closed = bool(constants.get("sourceInactiveTailConstantsClosed"))
    absorbable = bool(constants.get("absorbableByFiniteLowMidBlock"))
    inequality_closed = epsilon >= 0 and budget > 0 and slack > 0 and epsilon < budget
    theorem_closed = constants_closed and absorbable and inequality_closed

    data: dict[str, Any] = {
        "theoremName": "source-inactive tail constants consequence theorem",
        "proofClass": "symbolic identity",
        "basis": constants.get("basis"),
        "sourceWindow": constants.get("sourceWindow"),
        "globalTraceRatio": constants.get("globalTraceRatio"),
        "finiteSourceTop": constants.get("finiteSourceTop"),
        "finiteInactiveLambda3": constants.get("finiteInactiveLambda3"),
        "operatorPerturbationEta": constants.get("operatorPerturbationEta"),
        "absoluteTailBound": constants.get("absoluteTailBound"),
        "sourceTopLower": constants.get("sourceTopLower"),
        "normalizedEpsilonDelta": constants.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": constants.get("finiteLowMidSchurBudget"),
        "absorptionSlack": constants.get("absorptionSlack"),
        "epsilonOverBudget": constants.get("epsilonOverBudget"),
        "minMaxProofInCertifiedSourceModel": theorem_closed,
        "sourceInactiveTailConstantsClosed": theorem_closed,
        "sourceInactiveTailDominationConstantsClosed": theorem_closed,
        "absorbableByFiniteLowMidBlock": theorem_closed,
        "statuses": {
            "constantsInputStatus": status(
                "source-inactive constants input",
                constants_closed,
                "The detailed constants theorem supplies the normalized inactive-tail source constant.",
                blocker=None if constants_closed else "Close the detailed source-inactive constants theorem.",
            ),
            "absorptionInequalityStatus": status(
                "inactive tail absorbed by low/mid Schur budget",
                inequality_closed and absorbable,
                "The normalized inactive-tail constant is strictly below the finite low/mid Schur budget.",
                blocker=None
                if inequality_closed and absorbable
                else "Certify normalizedEpsilonDelta < finiteLowMidSchurBudget.",
            ),
            "tailConstantsConsequenceStatus": status(
                "terminal source-inactive tail constants consequence",
                theorem_closed,
                "The high-block theorem may use the exported constants and strict absorption slack.",
                blocker=None if theorem_closed else "Close the constants input or absorption inequality.",
            ),
        },
        "formalProof": [
            "Import the certified normalized inactive-tail constant.",
            "Import the certified finite low/mid Schur absorption budget.",
            "Check the strict scalar inequality normalizedEpsilonDelta < finiteLowMidSchurBudget.",
            "Export the constants as the terminal input contract for the high-block exhaustion theorem.",
        ],
        "nextProofTarget": (
            "Use this terminal constants consequence in high_block_exhaustion_theorem."
            if theorem_closed
            else "Close the failed source-inactive constants consequence status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Source-inactive tail constants consequence theorem")
    print(f"  epsilon: {epsilon:.12e}")
    print(f"  budget: {budget:.12e}")
    print(f"  slack: {slack:.12e}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
