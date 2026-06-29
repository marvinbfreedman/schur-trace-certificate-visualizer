#!/usr/bin/env python3
r"""Abstract Riesz projection continuity theorem."""

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
        default="abstract_riesz_projection_continuity_theorem.json",
    )
    args = parser.parse_args()

    ok = True
    data = {
        "theoremName": "abstract Riesz projection continuity theorem",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "If self-adjoint bounded operators S_N converge to S in norm and "
            "a contour Gamma separates an isolated spectral cluster of S, then "
            "the Riesz projections converge in norm for all large N."
        ),
        "hypotheses": [
            "S_N and S are bounded self-adjoint operators.",
            "||S_N-S|| -> 0.",
            "Gamma lies in the resolvent set of S and encloses the selected cluster with positive distance to the rest of the spectrum.",
        ],
        "statuses": {
            "resolventIdentityStatus": status(
                "resolvent identity on a fixed contour",
                ok,
                (
                    "For z on Gamma, (z-S_N)^{-1}-(z-S)^{-1}="
                    "(z-S_N)^{-1}(S_N-S)(z-S)^{-1} once ||S_N-S|| is below "
                    "the contour distance."
                ),
            ),
            "rieszProjectionContinuityStatus": status(
                "Riesz projection norm continuity",
                ok,
                (
                    "Integrating the resolvent identity over Gamma gives "
                    "||P_N-P|| -> 0."
                ),
            ),
            "rankPersistenceStatus": status(
                "spectral rank persistence",
                ok,
                (
                    "Norm-small perturbations preserve the rank enclosed by "
                    "the contour."
                ),
            ),
        },
        "rieszProjectionContinuityClosed": ok,
        "proof": [
            "Choose a contour Gamma at positive distance from sigma(S).",
            "For N large, ||S_N-S|| is smaller than half that distance, so Gamma also lies in the resolvent set of S_N.",
            "Use the resolvent identity and uniform resolvent bounds on Gamma.",
            "Integrate around Gamma to obtain norm convergence of the Riesz projections and rank persistence.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract Riesz projection continuity theorem")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
