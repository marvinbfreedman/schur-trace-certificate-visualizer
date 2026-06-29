#!/usr/bin/env python3
r"""Abstract compact-compression norm convergence theorem."""

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
        "--projection-json",
        default="abstract_mosco_projection_consequence_theorem.json",
    )
    parser.add_argument(
        "--finite-rank-json",
        default="abstract_finite_rank_compression_convergence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_compact_compression_norm_convergence_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.projection_json)
    finite_rank = load(args.finite_rank_json)
    projection_ok = bool(projection.get("strongProjectionConvergenceClosed"))
    finite_rank_ok = bool(finite_rank.get("finiteRankCompressionConvergenceClosed"))
    ok = projection_ok and finite_rank_ok

    data = {
        "theoremName": "abstract compact compression norm convergence theorem",
        "proofClass": "abstract functional analytic lemma",
        "projectionJson": args.projection_json,
        "finiteRankJson": args.finite_rank_json,
        "statement": (
            "If P_N -> P strongly are orthogonal projections and S is compact, "
            "then ||P_N S P_N - P S P|| -> 0 in operator norm."
        ),
        "hypotheses": [
            "P_N and P are orthogonal projections on one Hilbert space.",
            "P_N -> P strongly.",
            "S is a compact bounded operator.",
        ],
        "statuses": {
            "projectionInputStatus": status(
                "strong projection convergence input",
                projection_ok,
                "The Mosco projection theorem supplies P_N -> P strongly.",
            ),
            "finiteRankApproximationStatus": status(
                "finite-rank compact approximation",
                finite_rank_ok,
                (
                    "The abstract finite-rank compression theorem handles the "
                    "finite-dimensional approximation step."
                ),
            ),
            "compactCompressionNormConvergenceStatus": status(
                "compact compression norm convergence",
                ok,
                (
                    "The finite-rank approximation plus uniform boundedness of "
                    "projections gives ||P_NSP_N-PSP|| -> 0."
                ),
            ),
        },
        "compactCompressionNormConvergenceClosed": ok,
        "compactSourceNormConvergenceClosed": ok,
        "proof": [
            "Write P_NSP_N-PSP as terms involving (P_N-P)S and S(P_N-P).",
            "Approximate S by finite-rank F.",
            "Strong convergence of P_N is uniform on finite-dimensional subspaces, so the compressed finite-rank error tends to zero.",
            "The remaining terms are bounded by the compact approximation error.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Strong projection convergence or finite-rank compression input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract compact compression norm convergence theorem")
    print(f"  projection input: {projection_ok}")
    print(f"  finite-rank input: {finite_rank_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
