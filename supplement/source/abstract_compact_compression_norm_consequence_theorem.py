#!/usr/bin/env python3
r"""Terminal consequence for abstract compact-compression norm convergence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compression-json",
        default="abstract_compact_compression_norm_convergence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_compact_compression_norm_consequence_theorem.json",
    )
    args = parser.parse_args()

    compression = load(args.compression_json)
    theorem_closed = bool(
        compression.get("compactCompressionNormConvergenceClosed")
        and compression.get("compactSourceNormConvergenceClosed")
    )

    data: dict[str, Any] = {
        "theoremName": "abstract compact compression norm consequence theorem",
        "proofClass": "symbolic identity",
        "statement": (
            "The compact-compression norm theorem supplies the norm convergence "
            "input ||P_N S P_N - P S P|| -> 0 needed by downstream Riesz "
            "spectral projection arguments."
        ),
        "compactCompressionNormConvergenceClosed": theorem_closed,
        "compactSourceNormConvergenceClosed": theorem_closed,
        "statuses": {
            "compactCompressionInputStatus": status(
                "abstract compact-compression input",
                theorem_closed,
                "The detailed abstract compact-compression theorem has closed norm convergence.",
                blocker=None if theorem_closed else "Close abstract_compact_compression_norm_convergence_theorem.",
            ),
            "compactCompressionNormConsequenceStatus": status(
                "compact-compression norm convergence consequence",
                theorem_closed,
                "Downstream theorems may use the compact-source norm convergence flag.",
                blocker=None if theorem_closed else "Close the compact-compression input.",
            ),
        },
        "formalProof": [
            "Import the detailed abstract compact-compression norm theorem.",
            "Read off compactCompressionNormConvergenceClosed and compactSourceNormConvergenceClosed.",
            "Expose only the terminal norm-convergence flag to spectral-projection consumers.",
        ],
        "nextProofTarget": (
            "Use this terminal consequence in compact spectral projection theorems."
            if theorem_closed
            else "Close the compact-compression consequence status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Abstract compact compression norm consequence theorem")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
