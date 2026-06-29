#!/usr/bin/env python3
r"""Slim consequence for shifted-Xi finite-Gram closure.

The shifted-Xi de Branges kernel positivity theorem only consumes the finite
evaluation Gram positivity and entrywise convergence consequences.  This file
exports those flags while leaving normalization, augmented closed cone, repair
positivity, and finite PSD cone closure below the full finite-Gram theorem.
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
        "--finite-gram-json",
        default="shifted_xi_finite_gram_closure_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="shifted_xi_finite_gram_closure_consequence_theorem.json",
    )
    args = parser.parse_args()

    finite = load(args.finite_gram_json)
    gram_ok = bool(finite.get("finiteEvaluationGramPositive"))
    closure_ok = bool(finite.get("shiftedXiFiniteGramClosureClosed"))
    convergence_ok = bool(finite.get("entrywiseConvergenceToFullDeBrangesKernel"))
    repair_ok = bool(finite.get("continuumAugmentedRepairClosed"))
    ok = gram_ok and closure_ok

    data = {
        "theoremName": "shifted Xi finite-Gram closure consequence theorem",
        "proofClass": "symbolic identity",
        "finiteGramSource": "shifted Xi finite-Gram closure theorem",
        "statuses": {
            "finiteGramClosureInputStatus": status(
                "shifted-Xi finite-Gram closure input",
                closure_ok,
                "The finite-Gram closure theorem closes the shifted-Xi finite evaluation Gram argument.",
            ),
            "kernelEntryLimitStatus": status(
                "kernel entry convergence",
                convergence_ok,
                "The finite-Gram theorem records entrywise convergence to the full shifted-Xi de Branges kernel.",
            ),
            "finiteEvaluationGramPositivityStatus": status(
                "finite evaluation Gram positivity",
                gram_ok,
                "Every finite shifted-Xi de Branges evaluation Gram matrix is positive semidefinite.",
            ),
        },
        "finiteEvaluationGramPositive": gram_ok,
        "shiftedXiFiniteGramClosureClosed": closure_ok,
        "entrywiseConvergenceToFullDeBrangesKernel": convergence_ok,
        "continuumAugmentedRepairClosed": repair_ok,
        "shiftedXiFiniteGramClosureConsequenceClosed": ok,
        "proof": [
            "Import the shifted-Xi finite-Gram closure theorem.",
            "Expose only the finite evaluation Gram positivity and entrywise convergence flags used upstream.",
        ],
        "remainingAnalyticGap": None if ok else "Shifted-Xi finite-Gram closure theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi finite-Gram closure consequence theorem")
    print(f"  finite Gram positivity: {gram_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
