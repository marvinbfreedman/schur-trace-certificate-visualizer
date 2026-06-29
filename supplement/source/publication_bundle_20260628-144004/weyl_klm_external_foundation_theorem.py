#!/usr/bin/env python3
r"""External foundation theorem for Weyl/KLM normalization and kernel algebra.

This theorem aggregates the convention-level links that used to be checked by
searching prose in the notes and draft:

1. Riemann kernel/Xi normalization;
2. hbar=1 KLM/Weyl equivalence;
3. Weyl symbol to coordinate-kernel transport;
4. full-line to half-line parity reduction.
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
    parser.add_argument("--riemann-json", default="riemann_kernel_normalization_consequence_theorem.json")
    parser.add_argument("--klm-weyl-json", default="klm_weyl_hbar1_equivalence_theorem.json")
    parser.add_argument("--symbol-kernel-json", default="weyl_symbol_kernel_transport_consequence_theorem.json")
    parser.add_argument("--parity-json", default="parity_halfline_reduction_theorem.json")
    parser.add_argument("--json-out", default="weyl_klm_external_foundation_theorem.json")
    args = parser.parse_args()

    riemann = load(args.riemann_json)
    klm_weyl = load(args.klm_weyl_json)
    symbol_kernel = load(args.symbol_kernel_json)
    parity = load(args.parity_json)

    riemann_ok = bool(riemann.get("riemannKernelNormalizationClosed"))
    klm_ok = bool(klm_weyl.get("klmWeylHbar1EquivalenceClosed"))
    symbol_ok = bool(symbol_kernel.get("weylSymbolKernelTransportClosed"))
    parity_ok = bool(parity.get("parityHalflineReductionClosed"))
    foundation_ok = riemann_ok and klm_ok and symbol_ok and parity_ok

    data = {
        "theoremName": "Weyl/KLM external foundation theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "riemannKernelNormalizationStatus": status(
                "Riemann kernel formula and harmless Xi normalization",
                riemann_ok,
                "Imported Riemann kernel normalization theorem.",
                blocker=None if riemann_ok else "Close riemann_kernel_normalization_theorem.py.",
            ),
            "klmWeylEquivalenceStatus": status(
                "KLM convention equals Weyl positivity for hbar=1",
                klm_ok,
                "Imported hbar=1 KLM/Weyl equivalence theorem.",
                blocker=None if klm_ok else "Close klm_weyl_hbar1_equivalence_theorem.py.",
            ),
            "symbolKernelTransportStatus": status(
                "phase-space symbol transported to coordinate Weyl kernel",
                symbol_ok,
                "Imported Weyl symbol coordinate-kernel transport theorem.",
                blocker=None if symbol_ok else "Close weyl_symbol_kernel_transport_theorem.py.",
            ),
            "parityHalflineReductionStatus": status(
                "full-line Weyl kernel reduced to even/odd half-line sectors",
                parity_ok,
                "Imported parity half-line reduction theorem.",
                blocker=None if parity_ok else "Close parity_halfline_reduction_theorem.py.",
            ),
            "externalFoundationStatus": status(
                "external Weyl/KLM foundation",
                foundation_ok,
                (
                    "All convention and kernel-algebra links from the "
                    "normalized Volterra/Weyl certificate to the KLM-facing "
                    "Weyl positivity statement are theorem-backed."
                ),
                blocker=None if foundation_ok else "Close the open foundation theorem input.",
            ),
        },
        "externalFoundationClosed": foundation_ok,
        "klmWeylHbar1EquivalenceClosed": klm_ok,
        "proof": [
            "Import the Riemann kernel normalization theorem.",
            "Import the hbar=1 KLM/Weyl equivalence theorem.",
            "Import the Weyl symbol to coordinate-kernel transport theorem.",
            "Import the parity half-line reduction theorem.",
            "Combine them as the external normalization/kernel-algebra foundation.",
        ],
        "remainingAnalyticGap": None if foundation_ok else "One foundation theorem input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Weyl/KLM external foundation theorem")
    print(f"  Riemann normalization: {riemann_ok}")
    print(f"  KLM/Weyl hbar=1: {klm_ok}")
    print(f"  symbol-kernel transport: {symbol_ok}")
    print(f"  parity reduction: {parity_ok}")
    print(f"  foundation closed: {foundation_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
