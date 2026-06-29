#!/usr/bin/env python3
r"""Finite-Gram closure theorem for shifted-Xi de Branges kernels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--normalization-json",
        default="shifted_xi_debranges_normalization_theorem.json",
    )
    parser.add_argument(
        "--closed-cone-json",
        default="xi_augmented_closed_cone_limit_theorem.json",
    )
    parser.add_argument("--psd-closure-json", default="finite_psd_cone_closure_theorem.json")
    parser.add_argument("--json-out", default="shifted_xi_finite_gram_closure_theorem.json")
    args = parser.parse_args()

    normalization = load(args.normalization_json)
    cone = load(args.closed_cone_json)
    psd = load(args.psd_closure_json)

    normalization_ok = bool(normalization.get("shiftedXiDeBrangesNormalizationClosed"))
    finite_positive_ok = bool(cone.get("closedConeConclusionClosed"))
    entrywise_ok = bool(cone.get("entrywiseConvergenceToFullDeBrangesKernel"))
    repair_ok = bool(cone.get("continuumAugmentedRepairClosed"))
    psd_ok = bool(psd.get("finitePsdConeClosureClosed"))
    gram_ok = normalization_ok and finite_positive_ok and entrywise_ok and repair_ok and psd_ok

    data = {
        "theoremName": "shifted Xi finite-Gram closure theorem",
        "proofClass": "analytic proof",
        "normalizationJson": args.normalization_json,
        "closedConeJson": args.closed_cone_json,
        "psdClosureJson": args.psd_closure_json,
        "statement": (
            "For each finite upper-half-plane evaluation set and each "
            "0<omega<1/2, the shifted-Xi de Branges Gram matrix is the "
            "entrywise limit of positive finite augmented pullback Gram matrices."
        ),
        "statuses": {
            "normalizationInputStatus": status(
                "shifted-Xi normalization input",
                normalization_ok,
                "The target entries are exactly the shifted-Xi de Branges kernel entries.",
            ),
            "finitePositiveConeInputStatus": status(
                "finite positive augmented pullback input",
                finite_positive_ok,
                "The augmented closed-cone theorem supplies positive finite pullback Gram matrices.",
            ),
            "entrywiseKernelConvergenceInputStatus": status(
                "entrywise kernel convergence input",
                entrywise_ok,
                "K_{E,N}(w,z) converges entrywise to K_E(w,z) on finite evaluation sets.",
            ),
            "repairPositivityInputStatus": status(
                "continuum repair positivity input",
                repair_ok,
                "The augmented repair remains nonnegative in the transported trace space.",
            ),
            "finitePsdClosureInputStatus": status(
                "finite PSD cone closure input",
                psd_ok,
                "Entrywise limits of finite PSD matrices remain PSD.",
            ),
            "finiteEvaluationGramPositivityStatus": status(
                "finite evaluation Gram positivity",
                gram_ok,
                (
                    "For each finite evaluation set, the Gram matrix of "
                    "K_{E_omega} is an entrywise limit of PSD matrices and is therefore PSD."
                ),
            ),
        },
        "finiteEvaluationGramPositive": gram_ok,
        "shiftedXiFiniteGramClosureClosed": gram_ok,
        "entrywiseConvergenceToFullDeBrangesKernel": entrywise_ok,
        "continuumAugmentedRepairClosed": repair_ok,
        "proof": [
            "Fix omega and a finite set of points in C_+.",
            "Use the normalization theorem to identify the target matrix as the shifted-Xi de Branges Gram.",
            "Use the augmented closed-cone theorem to get positive finite truncation matrices.",
            "Use entrywise convergence of those matrices to the target finite Gram matrix.",
            "Apply finite PSD cone closure.",
        ],
        "remainingAnalyticGap": None if gram_ok else "A finite-Gram closure input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi finite-Gram closure theorem")
    print(f"  normalization: {normalization_ok}")
    print(f"  finite positive input: {finite_positive_ok}")
    print(f"  entrywise convergence: {entrywise_ok}")
    print(f"  finite Gram positivity: {gram_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
