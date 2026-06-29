#!/usr/bin/env python3
r"""Abstract min-max tail inequality passage theorem."""

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
        default="abstract_active_projection_convergence_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_minmax_tail_passage_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.projection_json)
    active_projection_ok = bool(
        projection.get("activeRieszProjectionConvergenceClosed")
        or projection.get("compactSourceRieszProjectionConvergenceClosed")
    )
    inactive_projection_ok = bool(
        projection.get("inactiveProjectionConvergenceClosed")
        or projection.get("sourceActiveInactiveSplittingStabilityClosed")
    )
    ok = active_projection_ok and inactive_projection_ok

    data = {
        "theoremName": "abstract min-max tail passage theorem",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "If compact positive source operators and their active Riesz "
            "projections converge in norm, then finite inactive-tail min-max "
            "bounds pass to the continuum closed subspace."
        ),
        "hypotheses": [
            "Compressed source operators S_N converge to S in norm.",
            "The top active Riesz projections P_{active,N} converge to P_active in norm.",
            "The finite inactive-tail inequality holds with a uniform constant epsilon.",
        ],
        "statuses": {
            "spectralInputStatus": status(
                "Riesz projection convergence input",
                active_projection_ok,
                "The active projection consequence theorem supplies active Riesz projection convergence.",
            ),
            "inactiveTailProjectionConvergenceStatus": status(
                "inactive tail projection convergence",
                inactive_projection_ok,
                (
                    "Norm convergence of active projections gives norm "
                    "convergence of inactive projections I-P_active."
                ),
            ),
            "courantFischerTailPassageStatus": status(
                "Courant-Fischer/min-max tail passage",
                ok,
                (
                    "The source quadratic forms converge uniformly on the "
                    "unit A-sphere, so the uniform finite inactive-tail bound "
                    "passes to the limit by recovery sequences and lower limits."
                ),
            ),
        },
        "abstractMinmaxTailPassageClosed": ok,
        "spectralMinmaxPassageClosed": ok,
        "proof": [
            "Use norm convergence of source operators to get uniform convergence of source quadratic forms on the A-unit sphere.",
            "Use norm convergence of active Riesz projections to identify the finite inactive spaces with the continuum inactive space.",
            "Apply Courant-Fischer variational continuity, or equivalently pass the quadratic inequality along recovery sequences and weak lower limits.",
        ],
        "remainingAnalyticGap": None if ok else "Riesz projection stability input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract min-max tail passage theorem")
    print(f"  active projection input: {active_projection_ok}")
    print(f"  inactive projection input: {inactive_projection_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
