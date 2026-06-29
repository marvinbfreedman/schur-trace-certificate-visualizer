#!/usr/bin/env python3
r"""Abstract quotient-norm liminf principle.

This is a model-free Hilbert-space lemma.  It does not assert any Volterra,
theta, Mellin, trace-convergence, or compactness input.

Lemma.  Let X be the transported quotient of a Hilbert graph space V by a
closed trace map R, with

    ||x||_X = inf { ||f||_V : R f = x }.

Suppose x_N are finite trace coordinates with near-minimizing representatives
f_N, and along a subsequence:

1. f_N converges weakly in V to f;
2. lower semicontinuity gives ||f||_V <= liminf ||f_N||_V;
3. trace convergence identifies R f = x.

Then

    ||x||_X <= liminf ||x_N||_{X_N}.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


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
        "--json-out",
        default="abstract_trace_quotient_liminf_principle_theorem.json",
    )
    args = parser.parse_args()

    ok = True
    data = {
        "theoremName": "abstract trace quotient liminf principle",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "Weak representative compactness, trace identification, and graph "
            "norm lower semicontinuity imply the quotient-norm lower bound."
        ),
        "statuses": {
            "abstractQuotientLiminfPrincipleStatus": status(
                "abstract quotient liminf principle",
                ok,
                (
                    "For any near-minimizing representatives f_N, weak graph "
                    "convergence and trace identification give an admissible "
                    "limit representative f for x.  The quotient norm of x is "
                    "bounded by ||f||_V, and lower semicontinuity bounds "
                    "||f||_V by liminf ||f_N||_V."
                ),
            ),
            "noModelSpecificInputStatus": status(
                "no model-specific input",
                ok,
                (
                    "The lemma is only the Hilbert quotient inequality.  All "
                    "compactness and trace-convergence hypotheses are supplied "
                    "by separate augmented trace theorems."
                ),
            ),
        },
        "abstractTraceQuotientLiminfPrincipleClosed": ok,
        "noVolterraOrThetaClaimed": True,
        "proof": [
            "Choose representatives f_N with ||f_N||_V within epsilon_N of the finite quotient norm.",
            "Use the assumed weak representative limit f and trace identification Rf=x.",
            "By the quotient definition, ||x||_X <= ||f||_V.",
            "By weak lower semicontinuity, ||f||_V <= liminf ||f_N||_V.",
            "Let epsilon_N tend to zero.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract trace quotient liminf principle")
    print(f"  principle closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
