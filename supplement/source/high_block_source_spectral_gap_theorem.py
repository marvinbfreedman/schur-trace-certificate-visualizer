#!/usr/bin/env python3
r"""High-block source spectral-gap theorem."""

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


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    if not item and isinstance(data.get("statuses"), dict):
        item = data["statuses"].get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-rank-json",
        default="full_theta_source_noncollapse_consequence_theorem.json",
    )
    parser.add_argument(
        "--tail-constants-json",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_source_spectral_gap_theorem.json",
    )
    args = parser.parse_args()

    source_rank = load(args.source_rank_json)
    tail = load(args.tail_constants_json)

    rank_ok = bool(source_rank.get("fullThetaSourceNoncollapseIntervalTheoremClosed"))
    projector_ok = closed(source_rank, "rieszProjectorStabilityStatus")
    gap = source_rank.get("spectralGapToComplementS8")
    lower_after_tail = source_rank.get("lowerBoundAfterContinuumAndTail")
    source_top_lower = tail.get("sourceTopLower")
    gap_positive = bool(gap is not None and gap > 0 and lower_after_tail is not None and lower_after_tail > 0)
    ok = rank_ok and projector_ok and gap_positive

    data = {
        "theoremName": "high-block source spectral-gap theorem",
        "proofClass": "interval/ball certificate",
        "sourceRankJson": args.source_rank_json,
        "tailConstantsJson": args.tail_constants_json,
        "statement": (
            "The full-theta source operator has a certified isolated top-two "
            "active spectral cluster separated from the inactive complement."
        ),
        "spectralGapToComplement": gap,
        "lowerBoundAfterContinuumAndTail": lower_after_tail,
        "sourceTopLower": source_top_lower,
        "statuses": {
            "sourceRankInputStatus": status(
                "full-theta source noncollapse interval theorem",
                rank_ok,
                "The interval theorem proves the active two-dimensional source rank persists.",
            ),
            "projectorGapInputStatus": status(
                "source Riesz projector gap input",
                projector_ok and gap_positive,
                (
                    "The source noncollapse theorem records a positive gap to "
                    "the complement and a positive lower bound after continuum "
                    "and theta-tail errors."
                ),
            ),
            "highBlockSourceSpectralGapStatus": status(
                "high-block source spectral gap",
                ok,
                "The certified source gap supplies the contour separation required by the abstract Riesz theorem.",
            ),
        },
        "highBlockSourceSpectralGapClosed": ok,
        "sourceSpectralGapClosed": ok,
        "proof": [
            "Import the full-theta source noncollapse interval theorem.",
            "Use its certified positive spectral gap to the complement and positive post-tail lower bound.",
            "Declare the active contour gap for the top-two source cluster in the completed high-block source model.",
        ],
        "remainingAnalyticGap": None if ok else "Source noncollapse interval theorem or positive source gap is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block source spectral-gap theorem")
    print(f"  source rank input: {rank_ok}")
    print(f"  projector input: {projector_ok}")
    print(f"  positive gap: {gap_positive}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
