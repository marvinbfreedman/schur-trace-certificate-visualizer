#!/usr/bin/env python3
r"""Abstract bounded-form descent to a quotient Hilbert space.

This is a model-free Hilbert-space lemma.  It does not assert any theta,
Volterra, Mellin, or finite Schur certificate statement.

Let R:V -> X be a closed trace map and let X be the transported quotient
Hilbert space closure of Ran(R), with quotient norm

    ||x||_X = inf { ||f||_V : Rf=x }.

If a quadratic form d is compatible with the trace equivalence relation and
bounded in the quotient norm,

    |d(Rf,Rg)| <= C ||Rf||_X ||Rg||_X,

then d descends uniquely to a bounded closed quadratic form on X.
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
        default="abstract_bounded_form_quotient_descent_theorem.json",
    )
    args = parser.parse_args()

    ok = True
    data = {
        "theoremName": "abstract bounded form quotient descent theorem",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "A trace-compatible quadratic form that is bounded in the "
            "transported quotient norm descends uniquely to a bounded closed "
            "quadratic form on the quotient Hilbert space."
        ),
        "statuses": {
            "abstractQuotientDescentPrincipleStatus": status(
                "abstract quotient descent principle",
                ok,
                (
                    "Trace compatibility makes the value independent of the "
                    "chosen representative.  Quotient-norm boundedness makes "
                    "the descended form continuous on the quotient Hilbert "
                    "space, hence closed."
                ),
            ),
            "noModelSpecificInputStatus": status(
                "no model-specific input",
                ok,
                (
                    "The lemma only uses the quotient Hilbert structure and "
                    "boundedness.  Finite repair compatibility and bounds are "
                    "supplied by separate augmented theorems."
                ),
            ),
        },
        "abstractBoundedFormQuotientDescentClosed": ok,
        "noVolterraOrThetaClaimed": True,
        "proof": [
            "If Rf=Rf', then f-f' lies in the null fiber; trace compatibility gives the same form value.",
            "Therefore the form is well-defined on equivalence classes in Ran(R).",
            "The quotient-norm bound gives a continuous bilinear/quadratic form on Ran(R).",
            "Continuous extension to the Hilbert closure X is unique and closed.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract bounded form quotient descent theorem")
    print(f"  descent principle closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
