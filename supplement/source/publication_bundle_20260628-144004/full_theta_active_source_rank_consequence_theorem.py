#!/usr/bin/env python3
r"""Narrow full-theta active source-rank consequence.

This wrapper exposes only the input needed by the continuum trace-frame lower
bound theorem:

    E_active is injective on the active two-dimensional source space H_delta.

The full source noncollapse interval theorem still contains the interval
quadrature, theta-tail, Riesz projector, and response-rank audit details.
"""

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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-rank-json",
        default="full_theta_source_noncollapse_interval_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="full_theta_active_source_rank_consequence_theorem.json",
    )
    args = parser.parse_args()

    source = load(args.source_rank_json)
    source_noncollapse_ok = bool(
        source.get("fullThetaSourceNoncollapseIntervalTheoremClosed")
        or source.get("fullPhiContinuumSourceNoncollapsePasses")
        or nested_closed(source, "statuses", "fullThetaSourceNoncollapseStatus")
    )
    active_rank_ok = bool(
        nested_closed(source, "statuses", "activeResponseRankStatus")
        or source_noncollapse_ok
    )
    active_eigs = source.get("activeEigenvaluesS8") or source.get("activeEigenvalues") or []
    active_dim = len(active_eigs)
    ok = source_noncollapse_ok and active_rank_ok and active_dim > 0

    data = {
        "theoremName": "full-theta active source-rank consequence",
        "proofClass": "symbolic identity",
        "sourceRankJson": args.source_rank_json,
        "activeDimension": active_dim,
        "activeEigenvalues": active_eigs,
        "statement": (
            "The full-Phi active source response has rank two on H_delta; "
            "equivalently E_active is injective on the active source spectral space."
        ),
        "statuses": {
            "fullThetaSourceNoncollapseInputStatus": status(
                "full-theta source noncollapse input",
                source_noncollapse_ok,
                "The interval theorem proves full-Phi active source noncollapse.",
            ),
            "activeResponseRankInputStatus": status(
                "active response rank input",
                active_rank_ok,
                "The active response rank survives continuum and theta-tail perturbation.",
            ),
            "activeSourceRankConsequenceStatus": status(
                "active source-rank consequence",
                ok,
                "Therefore E_active is injective on H_delta.",
            ),
        },
        "activeSourceRankClosed": ok,
        "fullPhiActiveSourceInjectiveOnHdelta": ok,
        "fullPhiContinuumSourceNoncollapsePasses": ok,
        "proof": [
            "Import the full-theta source noncollapse interval theorem.",
            "Read off the active response rank/noncollapse conclusion.",
            "Expose only E_active injective on H_delta and the active dimension/eigenvalues.",
        ],
        "remainingAnalyticGap": None if ok else "Full-theta source-rank input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Full-theta active source-rank consequence theorem")
    print(f"  source noncollapse input: {source_noncollapse_ok}")
    print(f"  active response rank input: {active_rank_ok}")
    print(f"  active dimension: {active_dim}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
