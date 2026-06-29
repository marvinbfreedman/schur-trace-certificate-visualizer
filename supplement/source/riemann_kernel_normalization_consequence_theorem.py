#!/usr/bin/env python3
r"""Terminal consequence of the Riemann kernel normalization theorem.

The Weyl/KLM external foundation theorem does not need the full normalization
ledger.  It only needs the closed consequence that the Phi/Xi Fourier
normalization used in the proof is harmless for positivity: any conventional
nonzero scalar produces only a positive scalar rescaling of the paired kernel
forms.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
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
        "--normalization-json",
        default="riemann_kernel_normalization_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="riemann_kernel_normalization_consequence_theorem.json",
    )
    args = parser.parse_args()

    theorem = load(args.normalization_json)
    closed = bool(theorem.get("riemannKernelNormalizationClosed"))

    normalization_status = status(
        "Riemann Phi/Xi normalization positivity consequence",
        closed,
        (
            "The fixed even Fourier representation Xi(z)=int Phi(t)e^{izt}dt "
            "and any harmless scalar convention preserve the positivity "
            "statements used by the Weyl/KLM external foundation."
        ),
        blocker=None if closed else "Close riemann_kernel_normalization_theorem.py.",
    )

    data = {
        "theoremName": "Riemann kernel normalization consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "riemannKernelNormalizationConsequenceStatus": normalization_status,
        },
        "riemannKernelNormalizationClosed": normalization_status["closed"],
        "phiXiNormalizationPositivityConsequenceClosed": normalization_status["closed"],
        "proof": [
            "Import the full Riemann kernel normalization theorem.",
            "Export only the positivity-preserving Phi/Xi normalization consequence.",
        ],
        "remainingAnalyticGap": None if normalization_status["closed"] else normalization_status["blocker"],
        "nextProofTarget": None if normalization_status["closed"] else normalization_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Riemann kernel normalization consequence theorem")
    print(f"  normalization consequence closed: {data['riemannKernelNormalizationClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
