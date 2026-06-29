#!/usr/bin/env python3
r"""Narrow consequence for shifted-Xi de Branges kernel positivity.

The Hermite--Biehler endpoint passage only needs one exported fact:

    K_{E_omega} >= 0 on finite upper-half-plane evaluation sets,
    E_omega(z)=Xi(z+i omega), 0<omega<1/2.

This file reads the full shifted-Xi kernel positivity theorem and exposes only
that terminal positivity consequence.  Finite-Gram closure, entrywise limits,
and augmented trace repair remain below the full theorem.
"""

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
        "--kernel-positivity-json",
        default="shifted_xi_debranges_kernel_positivity_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="shifted_xi_debranges_kernel_positivity_consequence_theorem.json",
    )
    args = parser.parse_args()

    kernel = load(args.kernel_positivity_json)
    kernel_ok = bool(kernel.get("shiftedXiDeBrangesKernelPositiveClosed"))
    finite_gram_ok = bool(kernel.get("finiteEvaluationGramPositive"))
    repair_ok = bool(kernel.get("continuumAugmentedRepairClosed"))

    data = {
        "theoremName": "shifted Xi de Branges kernel positivity consequence theorem",
        "proofClass": "symbolic identity",
        "source": "shifted Xi de Branges kernel positivity theorem",
        "omegaRange": "0 < omega < 1/2",
        "kernel": "K_{E_omega}, E_omega(z)=Xi(z+i omega)",
        "statuses": {
            "shiftedXiKernelPositivityInputStatus": status(
                "shifted-Xi de Branges kernel positivity input",
                kernel_ok,
                (
                    "The full shifted-Xi kernel positivity theorem proves "
                    "positive definiteness on finite upper-half-plane "
                    "evaluation sets."
                ),
            ),
            "finiteGramSummaryStatus": status(
                "finite Gram positivity summary",
                finite_gram_ok,
                "The lower theorem records the finite evaluation Gram positivity input.",
            ),
            "augmentedRepairSummaryStatus": status(
                "augmented repair closure summary",
                repair_ok,
                "The lower theorem records the continuum augmented repair closure input.",
            ),
        },
        "shiftedXiDeBrangesKernelPositiveConsequenceClosed": kernel_ok,
        "shiftedXiDeBrangesKernelPositiveClosed": kernel_ok,
        "finiteEvaluationGramPositive": finite_gram_ok,
        "continuumAugmentedRepairClosed": repair_ok,
        "proof": [
            "Import the full shifted-Xi de Branges kernel positivity theorem.",
            "Export only the terminal positivity statement used by endpoint zero exclusion.",
        ],
        "notExportedHere": [
            "finite theta pullback construction",
            "augmented trace convergence",
            "repair positivity",
            "closed-cone limiting argument",
        ],
        "remainingAnalyticGap": None
        if kernel_ok
        else "Shifted-Xi de Branges kernel positivity theorem is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi de Branges kernel positivity consequence theorem")
    print(f"  kernel positivity consequence: {kernel_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
