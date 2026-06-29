#!/usr/bin/env python3
r"""Positive Volterra tail form theorem.

This theorem isolates the exact input needed by the augmented repair transport
layer.  The repair transport does not need the full original Weyl/KLM bridge;
it only needs nonnegativity of the diagonal Volterra tail form on the completed
Mellin/Volterra trace-fiber domain.

Let

    P(f)=||G_+ f||^2,       M(f)=||G_- f||^2.

The closed Green-lift theorem gives the completed-domain branch contraction

    G_- = C K E G_+,        ||C K E|| <= 1,

on the Green-minimizer trace image.  Therefore

    M(f) <= P(f),
    T_volt(f)=P(f)-M(f) >= 0.

The statement is independent of omega and sits below the original coordinate
Weyl/KLM packaging.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--green-closure-json",
        default="continuum_green_lift_closure_theorem.json",
    )
    parser.add_argument("--json-out", default="volterra_tail_positive_form_theorem.json")
    args = parser.parse_args()

    green = load(args.green_closure_json)
    green_ok = bool(
        green.get("continuumGreenLiftClosureClosed")
        or green.get("closedOnCompletedVolterraDomain")
        or nested_closed(green, "statuses", "greenLiftContractionStatus")
    )
    tail_ok = green_ok

    data = {
        "theoremName": "positive Volterra tail form theorem",
        "proofClass": "analytic proof",
        "greenClosureJson": args.green_closure_json,
        "statuses": {
            "completedGreenLiftContractionStatus": status(
                "completed Green-lift contraction",
                green_ok,
                (
                    "The imported continuum Green-lift closure theorem proves "
                    "the compressed branch contraction ||C K E||<=1 on the "
                    "completed trace-fiber domain."
                ),
                blocker=None if green_ok else "Close continuum_green_lift_closure_theorem.py.",
            ),
            "volterraTailPositivityStatus": status(
                "positive Volterra tail form",
                tail_ok,
                (
                    "Since G_- = C K E G_+ and ||C K E||<=1, one has "
                    "||G_-f||^2<=||G_+f||^2.  Thus the diagonal Volterra tail "
                    "form P-M is nonnegative on the completed domain."
                ),
                blocker=None if tail_ok else "Need the completed Green-lift contraction.",
            ),
        },
        "volterraTailPositiveFormClosed": tail_ok,
        "omegaIndependent": True,
        "notClaimedHere": (
            "This theorem does not assert original coordinate Weyl positivity "
            "or the KLM positive-type condition.  Those are packaged by the "
            "separate uniform omega Weyl/KLM bridge."
        ),
        "proof": [
            "Import the completed-domain Green-lift contraction ||C K E||<=1.",
            "Use the branch identity G_-=C K E G_+ on the Green-minimizer trace image.",
            "Conclude ||G_-f||^2<=||G_+f||^2.",
            "Therefore the diagonal Volterra tail form P-M is nonnegative.",
        ],
        "remainingAnalyticGap": None if tail_ok else "The Green-lift contraction input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Positive Volterra tail form theorem")
    print(f"  Green-lift contraction: {green_ok}")
    print(f"  tail positivity: {tail_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
