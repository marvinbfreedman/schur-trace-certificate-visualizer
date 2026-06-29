#!/usr/bin/env python3
r"""Abstract closed-LSC nonnegative cone principle.

This is a pure functional-analytic lemma used by the augmented repair limit.
It is deliberately independent of the Volterra, trace, quotient, and finite
certificate machinery.

Lemma.  Let q be the closed lower-semicontinuous envelope of a family of
quadratic forms q_N on a recovery core.  If q_N(f_N) >= 0 for every core vector
in every approximating core, then q(f) >= 0 for every f in the completed form
domain.

Equivalently, the positive cone of closed lower-semicontinuous quadratic forms
is closed under Mosco/lower-envelope limits.
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
        default="closed_lsc_nonnegative_cone_principle_theorem.json",
    )
    args = parser.parse_args()

    principle_ok = True
    data = {
        "theoremName": "closed LSC nonnegative cone principle",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "If a closed lower-semicontinuous quadratic form is the Mosco/"
            "lower-envelope limit of nonnegative core quadratic forms, then "
            "the completed closed form is nonnegative."
        ),
        "statuses": {
            "closedLscConePrincipleStatus": status(
                "closed LSC nonnegative cone principle",
                principle_ok,
                (
                    "For any completed-domain vector f, the closed lower "
                    "envelope is computed from neighborhood liminf infima of "
                    "the core forms.  Since every core form is nonnegative, "
                    "each such neighborhood infimum is nonnegative, so the "
                    "closed envelope cannot assign a negative value."
                ),
            ),
            "noModelSpecificInputStatus": status(
                "no model-specific input",
                principle_ok,
                (
                    "This lemma uses only closed lower-semicontinuity and the "
                    "ordered cone structure of quadratic forms; all Volterra/"
                    "trace assumptions are supplied by separate theorems."
                ),
            ),
        },
        "closedLscNonnegativeConePrincipleClosed": principle_ok,
        "noTraceOrQuotientClaimed": True,
        "proof": [
            "The pointwise nonnegative core forms lie in the closed convex cone of nonnegative quadratic forms.",
            "The closed lower envelope is the Gamma/Mosco liminf envelope: a supremum of neighborhood liminf infima.",
            "Every neighborhood liminf infimum of nonnegative core forms is nonnegative.",
            "A negative value of the closed envelope would therefore require a negative core value in some neighborhood, contradicting core nonnegativity.",
            "Therefore the lower-semicontinuous completed envelope remains in the nonnegative cone.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed LSC nonnegative cone principle")
    print(f"  principle closed: {principle_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
