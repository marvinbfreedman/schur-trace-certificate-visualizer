#!/usr/bin/env python3
r"""Riemann kernel normalization theorem.

This theorem isolates the harmless normalization link between the Riemann
kernel Phi and the shifted Xi transform used by the Weyl/KLM ledger.  A
nonzero scalar multiple in the Fourier representation of Xi rescales the
resulting Gram kernels by a positive scalar in the positivity statements, so
it cannot change positive semidefiniteness.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", default="riemann_kernel_normalization_theorem.json")
    args = parser.parse_args()

    data = {
        "theoremName": "Riemann kernel normalization theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "riemannKernelFormulaStatus": status(
                "Riemann kernel Fourier representation",
                True,
                (
                    "The even Riemann kernel Phi is the Fourier kernel for Xi "
                    "in the normalization used throughout the Weyl/KLM ledger."
                ),
            ),
            "harmlessScalarNormalizationStatus": status(
                "harmless Xi scalar normalization",
                True,
                (
                    "Multiplying Xi or Phi by a fixed nonzero scalar only "
                    "rescales the finite Gram matrices by a positive scalar "
                    "after the paired conjugation used in the kernel."
                ),
            ),
        },
        "riemannKernelNormalizationClosed": True,
        "proof": [
            "Use the fixed even Fourier representation Xi(z)=int Phi(t)e^{izt}dt.",
            "Track any conventional scalar as a global kernel rescaling.",
            "Global positive scalar rescalings preserve positive semidefiniteness.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Riemann kernel normalization theorem")
    print("  normalization closed: True")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
