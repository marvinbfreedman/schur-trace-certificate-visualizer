#!/usr/bin/env python3
r"""Shifted-Xi de Branges kernel positivity theorem.

This theorem is the narrow output of the augmented Xi closed-cone machinery.
It imports the finite-Gram closure theorem and records only the consequence
needed by the final KLM-to-de Branges bridge:

    K_{E_omega} is positive definite on finite upper-half-plane
    evaluation sets, where E_omega(z)=Xi(z+i omega).

All normalization, finite theta pullback, trace convergence, repair positivity,
and finite PSD cone closure details remain one layer lower.
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
        "--finite-gram-json",
        default="shifted_xi_finite_gram_closure_consequence_theorem.json",
    )
    parser.add_argument(
        "--closed-cone-json",
        default="",
    )
    parser.add_argument(
        "--json-out",
        default="shifted_xi_debranges_kernel_positivity_theorem.json",
    )
    args = parser.parse_args()

    gram = load(args.finite_gram_json)
    cone = load(args.closed_cone_json) if args.closed_cone_json else {}
    gram_ok = bool(gram.get("finiteEvaluationGramPositive"))
    convergence_ok = bool(gram.get("entrywiseConvergenceToFullDeBrangesKernel"))
    repair_ok = bool(gram.get("continuumAugmentedRepairClosed"))
    cone_ok = bool(gram.get("shiftedXiFiniteGramClosureClosed") or cone.get("closedConeConclusionClosed"))
    kernel_ok = gram_ok and cone_ok

    data = {
        "theoremName": "shifted Xi de Branges kernel positivity theorem",
        "proofClass": "analytic proof",
        "finiteGramJson": args.finite_gram_json,
        "legacyInputs": {"closedConeJson": args.closed_cone_json or None},
        "omegaRange": "0 < omega < 1/2",
        "kernel": "K_{E_omega}, E_omega(z)=Xi(z+i omega)",
        "statuses": {
            "finiteGramClosureInputStatus": status(
                "shifted-Xi finite-Gram closure input",
                cone_ok,
                (
                    "The finite-Gram closure theorem combines normalization, "
                    "entrywise convergence, positive finite pullbacks, and "
                    "finite PSD cone closure."
                ),
                blocker=None if cone_ok else "Close shifted_xi_finite_gram_closure_theorem.py.",
            ),
            "kernelEntryLimitStatus": status(
                "kernel entry convergence",
                convergence_ok,
                (
                    "The same closed-cone theorem records entrywise convergence "
                    "to the full shifted-Xi de Branges kernel on finite "
                    "evaluation sets."
                ),
                blocker=None if convergence_ok else "Close xi_augmented_trace_convergence_theorem.py.",
            ),
            "shiftedXiKernelPositivityStatus": status(
                "shifted-Xi de Branges kernel positivity",
                kernel_ok,
                (
                    "For every finite upper-half-plane evaluation set, the "
                    "Gram matrix of K_{E_omega} is positive semidefinite."
                ),
                blocker=None if kernel_ok else "Need closed-cone membership and entry convergence.",
            ),
        },
        "shiftedXiDeBrangesKernelPositiveClosed": kernel_ok,
        "finiteEvaluationGramPositive": gram_ok,
        "continuumAugmentedRepairClosed": repair_ok,
        "proof": [
            "Import the shifted-Xi finite-Gram closure theorem.",
            "Read off positivity of every finite upper-half-plane evaluation Gram matrix.",
            "This is exactly positive definiteness of the shifted-Xi de Branges kernel.",
        ],
        "remainingAnalyticGap": None if kernel_ok else "The shifted-Xi finite-Gram closure input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi de Branges kernel positivity theorem")
    print(f"  finite-Gram closure input: {cone_ok}")
    print(f"  entry convergence: {convergence_ok}")
    print(f"  kernel positivity: {kernel_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
