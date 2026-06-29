#!/usr/bin/env python3
r"""High-block compact-compression input theorem."""

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
        default="high_block_mosco_projection_consequence_theorem.json",
    )
    parser.add_argument(
        "--source-compactness-json",
        default="high_block_source_operator_compactness_consequence_theorem.json",
    )
    parser.add_argument(
        "--compression-json",
        default="abstract_compact_compression_norm_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_compact_compression_input_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.projection_json)
    source_compactness = load(args.source_compactness_json)
    compression = load(args.compression_json)

    projection_ok = bool(
        projection.get("highBlockMoscoProjectionClosed")
        or projection.get("strongProjectionConvergenceClosed")
    )
    source_compact_ok = bool(source_compactness.get("compactSourceOperatorClosed"))
    compression_ok = bool(compression.get("compactCompressionNormConvergenceClosed"))
    ok = projection_ok and source_compact_ok and compression_ok

    data = {
        "theoremName": "high-block compact-compression input theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The high-block source operator satisfies the hypotheses of the "
            "abstract compact-compression theorem, hence the compressed source "
            "operators converge in norm."
        ),
        "statuses": {
            "projectionInputStatus": status(
                "projection convergence input",
                projection_ok,
                "Mosco projection convergence supplies P_N -> P strongly.",
            ),
            "sourceCompactnessInputStatus": status(
                "source compactness input",
                source_compact_ok,
                "The high-block source compactness theorem supplies compactness of S=E^*E.",
            ),
            "abstractCompressionInputStatus": status(
                "abstract compact-compression theorem input",
                compression_ok,
                "The abstract theorem gives ||P_NSP_N-PSP|| -> 0 once its hypotheses are met.",
            ),
            "highBlockCompactCompressionStatus": status(
                "high-block compact compression",
                ok,
                "The projection and compactness inputs match the abstract compression theorem.",
            ),
        },
        "highBlockCompactCompressionClosed": ok,
        "compactCompressionNormConvergenceClosed": ok,
        "compactSourceNormConvergenceClosed": ok,
        "proof": [
            "Identify the Hilbert space with the completed A-high-block space.",
            "Use Mosco projection convergence for the finite closed subspaces.",
            "Use the source-range Hardy/Green compactness theorem for S=E^*E.",
            "Apply the abstract compact-compression theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Projection, source compactness, or abstract compression input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block compact-compression input theorem")
    print(f"  projection input: {projection_ok}")
    print(f"  source compactness input: {source_compact_ok}")
    print(f"  abstract compression input: {compression_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
