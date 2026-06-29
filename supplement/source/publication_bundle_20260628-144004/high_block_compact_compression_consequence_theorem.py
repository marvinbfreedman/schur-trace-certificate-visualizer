#!/usr/bin/env python3
r"""Narrow high-block compact-compression consequence theorem.

This wrapper exposes only the compact-compression facts needed by the
high-block tail-estimate continuum-passage theorem:

    high-block compact compression,
    compact-compression norm convergence,
    compact source norm convergence.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compact-compression-json",
        default="high_block_compact_compression_input_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_compact_compression_consequence_theorem.json",
    )
    args = parser.parse_args()

    compression = load(args.compact_compression_json)
    compact_ok = bool(compression.get("highBlockCompactCompressionClosed"))
    norm_ok = bool(compression.get("compactCompressionNormConvergenceClosed"))
    source_norm_ok = bool(compression.get("compactSourceNormConvergenceClosed"))
    ok = compact_ok and norm_ok and source_norm_ok

    data = {
        "theoremName": "high-block compact-compression consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The high-block compact-compression input theorem supplies the "
            "compact source compression consequences needed by the continuum "
            "high-block tail estimate passage."
        ),
        "statuses": {
            "highBlockCompactCompressionStatus": status(
                "high-block compact compression",
                compact_ok,
                "Imported from the high-block compact-compression input theorem.",
            ),
            "compactCompressionNormConvergenceStatus": status(
                "compact-compression norm convergence",
                norm_ok,
                "Imported from the high-block compact-compression input theorem.",
            ),
            "compactSourceNormConvergenceStatus": status(
                "compact source norm convergence",
                source_norm_ok,
                "Imported from the high-block compact-compression input theorem.",
            ),
            "compactCompressionConsequenceStatus": status(
                "compact-compression consequence",
                ok,
                (
                    "Together these consequences give the compact source "
                    "compression input used by the high-block continuum tail "
                    "passage."
                ),
            ),
        },
        "highBlockCompactCompressionClosed": compact_ok,
        "compactCompressionNormConvergenceClosed": norm_ok,
        "compactSourceNormConvergenceClosed": source_norm_ok,
        "highBlockCompactCompressionConsequenceClosed": ok,
        "proof": [
            "Import the high-block compact-compression input theorem.",
            "Extract the three compact-compression consequences used upstream.",
            "Expose only these consequences to the high-block tail-passage theorem.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "High-block compact compression, norm convergence, or source norm convergence is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block compact-compression consequence theorem")
    print(f"  high-block compact compression: {compact_ok}")
    print(f"  compact-compression norm convergence: {norm_ok}")
    print(f"  compact source norm convergence: {source_norm_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
