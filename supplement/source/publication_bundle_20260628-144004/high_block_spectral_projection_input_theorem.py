#!/usr/bin/env python3
r"""High-block spectral projection input theorem."""

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
        "--compact-compression-json",
        default="high_block_compact_compression_input_theorem.json",
    )
    parser.add_argument(
        "--source-gap-json",
        default="high_block_source_spectral_gap_theorem.json",
    )
    parser.add_argument(
        "--abstract-spectral-json",
        default="abstract_compact_source_spectral_projection_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_spectral_projection_input_theorem.json",
    )
    args = parser.parse_args()

    compact = load(args.compact_compression_json)
    gap = load(args.source_gap_json)
    spectral = load(args.abstract_spectral_json)

    compact_ok = bool(compact.get("highBlockCompactCompressionClosed"))
    gap_ok = bool(gap.get("highBlockSourceSpectralGapClosed"))
    spectral_ok = bool(spectral.get("compactSourceRieszProjectionConvergenceClosed"))
    ok = compact_ok and gap_ok and spectral_ok

    data = {
        "theoremName": "high-block spectral projection input theorem",
        "proofClass": "symbolic identity",
        "compactCompressionJson": args.compact_compression_json,
        "sourceGapJson": args.source_gap_json,
        "abstractSpectralJson": args.abstract_spectral_json,
        "statement": (
            "The high-block source model satisfies the norm-convergence and "
            "spectral-gap hypotheses of the abstract Riesz projection theorem."
        ),
        "statuses": {
            "compactCompressionInputStatus": status(
                "compact-compression norm convergence input",
                compact_ok,
                "The high-block compact-compression theorem gives norm convergence of compressed source operators.",
            ),
            "sourceGapInputStatus": status(
                "source spectral-gap input",
                gap_ok,
                "The high-block source spectral-gap theorem gives an isolated active top-two cluster.",
            ),
            "abstractSpectralInputStatus": status(
                "abstract spectral projection input",
                spectral_ok,
                "The abstract Riesz theorem gives active/inactive spectral projection stability.",
            ),
            "highBlockSpectralProjectionStatus": status(
                "high-block spectral projection stability",
                ok,
                "The high-block inputs satisfy the abstract spectral projection theorem.",
            ),
        },
        "highBlockSpectralProjectionClosed": ok,
        "compactSourceRieszProjectionConvergenceClosed": ok,
        "sourceActiveInactiveSplittingStabilityClosed": ok,
        "proof": [
            "Use compact-compression norm convergence for the compressed source operators.",
            "Use the interval-certified source spectral gap for the top-two active cluster.",
            "Apply the abstract Riesz projection continuity theorem through the abstract compact-source spectral theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Compact-compression, source-gap, or abstract spectral input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block spectral projection input theorem")
    print(f"  compact-compression input: {compact_ok}")
    print(f"  source-gap input: {gap_ok}")
    print(f"  abstract spectral input: {spectral_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
