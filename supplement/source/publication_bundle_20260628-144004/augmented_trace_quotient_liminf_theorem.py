#!/usr/bin/env python3
r"""Quotient-norm liminf theorem for augmented traces."""

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
        "--principle-json",
        default="abstract_trace_quotient_liminf_principle_theorem.json",
    )
    parser.add_argument(
        "--representative-compactness-json",
        default="augmented_trace_liminf_representative_compactness_consequence_theorem.json",
    )
    parser.add_argument("--liminf-json", default="")
    parser.add_argument("--trace-convergence-json", default="")
    parser.add_argument("--json-out", default="augmented_trace_quotient_liminf_theorem.json")
    args = parser.parse_args()

    principle = load(args.principle_json)
    representative = load(args.representative_compactness_json)

    principle_ok = bool(principle.get("abstractTraceQuotientLiminfPrincipleClosed"))
    representative_ok = bool(
        representative.get("representativeCompactnessAndTraceIdentificationClosed")
    )
    quotient_ok = principle_ok and representative_ok

    data = {
        "theoremName": "augmented trace quotient liminf theorem",
        "proofClass": "analytic proof",
        "principleJson": args.principle_json,
        "representativeCompactnessJson": args.representative_compactness_json,
        "legacyLiminfJson": args.liminf_json or None,
        "legacyTraceConvergenceJson": args.trace_convergence_json or None,
        "statement": (
            "If x_N=R_aug,N f_N converges to x and the finite quotient norms are "
            "bounded, then x belongs to X_aug and ||x||_{X_aug} <= liminf ||x_N||_{X_aug,N}."
        ),
        "statuses": {
            "abstractQuotientLiminfPrincipleInputStatus": status(
                "abstract quotient liminf principle input",
                principle_ok,
                (
                    "The abstract Hilbert quotient lemma converts weak "
                    "representative compactness, trace identification, and "
                    "graph lower semicontinuity into the quotient lower bound."
                ),
            ),
            "representativeCompactnessInputStatus": status(
                "representative compactness input",
                representative_ok,
                (
                    "The augmented representative theorem supplies bounded "
                    "near-minimizer compactness and trace limit identification."
                ),
            ),
            "traceQuotientLiminfStatus": status(
                "trace quotient liminf",
                quotient_ok,
                (
                    "The representative theorem verifies the hypotheses of "
                    "the abstract quotient-liminf principle in the augmented "
                    "trace setting."
                ),
            ),
        },
        "traceQuotientLiminfClosed": quotient_ok,
        "finiteQuotientLowerBoundClosed": quotient_ok,
        "proof": [
            "Import the abstract quotient-norm liminf principle.",
            "Import augmented representative compactness and trace identification.",
            "Apply the abstract principle to the verified augmented hypotheses.",
        ],
        "remainingAnalyticGap": None
        if quotient_ok
        else "The abstract quotient principle or representative compactness input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient liminf theorem")
    print(f"  abstract principle: {principle_ok}")
    print(f"  representative compactness: {representative_ok}")
    print(f"  quotient liminf: {quotient_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
