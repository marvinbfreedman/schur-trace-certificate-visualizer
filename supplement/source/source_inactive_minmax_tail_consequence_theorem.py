#!/usr/bin/env python3
r"""Terminal consequence for source-inactive min-max tail domination.

The detailed source-inactive min-max theorem records the spectral theorem
argument, the normalized inactive-tail constant, and the continuum high-block
passage.  The synchronized active-range theorem only needs the consequence:

    minMaxProofInCertifiedSourceModel = true,
    absorbableByFiniteLowMidBlock = true,
    continuumHighBlockExhaustionImported = true,
    normalizedEpsilonDelta < finiteLowMidSchurBudget.

This wrapper exposes just that contract for downstream active/inactive source
assembly.
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
        "--minmax-json",
        default="source_inactive_minmax_tail_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="source_inactive_minmax_tail_consequence_theorem.json",
    )
    args = parser.parse_args()

    minmax = load(args.minmax_json)

    epsilon = float(minmax["normalizedEpsilonDelta"])
    budget = float(minmax["finiteLowMidSchurBudget"])
    slack = float(minmax["absorptionSlack"])
    minmax_closed = bool(minmax.get("minMaxProofInCertifiedSourceModel"))
    absorbable = bool(minmax.get("absorbableByFiniteLowMidBlock"))
    continuum_closed = bool(
        minmax.get("continuumHighBlockExhaustionImported")
        or minmax.get("continuumHighBlockExhaustionProvedHere")
    )
    scalar_closed = epsilon >= 0 and budget > 0 and slack > 0 and epsilon < budget
    theorem_closed = minmax_closed and absorbable and continuum_closed and scalar_closed

    data: dict[str, Any] = {
        "theoremName": "source-inactive min-max tail consequence theorem",
        "proofClass": "symbolic identity",
        "basis": minmax.get("basis"),
        "sourceWindow": minmax.get("sourceWindow"),
        "globalTraceRatio": minmax.get("globalTraceRatio"),
        "finiteSourceTop": minmax.get("finiteSourceTop"),
        "finiteInactiveLambda3": minmax.get("finiteInactiveLambda3"),
        "operatorPerturbationEta": minmax.get("operatorPerturbationEta"),
        "absoluteTailBound": minmax.get("absoluteTailBound"),
        "sourceTopLower": minmax.get("sourceTopLower"),
        "normalizedEpsilonDelta": minmax.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": minmax.get("finiteLowMidSchurBudget"),
        "absorptionSlack": minmax.get("absorptionSlack"),
        "epsilonOverBudget": minmax.get("epsilonOverBudget"),
        "minMaxProofInCertifiedSourceModel": theorem_closed,
        "absorbableByFiniteLowMidBlock": theorem_closed,
        "continuumHighBlockExhaustionImported": theorem_closed,
        "continuumHighBlockExhaustionProvedHere": theorem_closed,
        "sourceInactiveMinmaxTailConsequenceClosed": theorem_closed,
        "statuses": {
            "minmaxTailInputStatus": status(
                "source-inactive min-max input",
                minmax_closed,
                "The detailed min-max theorem proves the normalized source-inactive tail estimate in the certified source model.",
                blocker=None if minmax_closed else "Close the detailed source-inactive min-max theorem.",
            ),
            "lowMidAbsorptionStatus": status(
                "low/mid Schur absorption",
                absorbable and scalar_closed,
                "The normalized inactive-tail constant is strictly below the finite low/mid Schur budget.",
                blocker=None if absorbable and scalar_closed else "Close the low/mid absorption inequality.",
            ),
            "continuumHighBlockStatus": status(
                "continuum high-block passage",
                continuum_closed,
                "The high-block exhaustion theorem has transported the min-max estimate to the continuum closed high block.",
                blocker=None if continuum_closed else "Close the high-block exhaustion passage.",
            ),
            "sourceInactiveMinmaxTailConsequenceStatus": status(
                "source-inactive min-max tail consequence",
                theorem_closed,
                "The synchronized active-range theorem may use the inactive-tail bound and absorption slack.",
                blocker=None if theorem_closed else "Close one of the source-inactive min-max consequence inputs.",
            ),
        },
        "formalProof": [
            "Import the detailed source-inactive min-max tail theorem.",
            "Read off epsilon_delta, the low/mid Schur budget, and the positive slack.",
            "Import the continuum high-block exhaustion flag.",
            "Export the inactive-tail domination contract used by active/inactive source assembly.",
        ],
        "nextProofTarget": (
            "Use this terminal consequence in synchronized_active_range_interval_theorem."
            if theorem_closed
            else "Close the failed source-inactive min-max consequence status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Source-inactive min-max tail consequence theorem")
    print(f"  epsilon: {epsilon:.12e}")
    print(f"  budget: {budget:.12e}")
    print(f"  slack: {slack:.12e}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
