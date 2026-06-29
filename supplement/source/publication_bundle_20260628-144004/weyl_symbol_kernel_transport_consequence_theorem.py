#!/usr/bin/env python3
r"""Terminal consequence of the Weyl symbol coordinate-kernel transport theorem.

The external foundation theorem only needs the endpoint statement that, in the
fixed Weyl convention, positivity of the phase-space Weyl operator
Op^W(sigma_omega) is equivalently positivity of the associated coordinate
kernel quadratic form on the Weyl test domain.
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
        "--transport-json",
        default="weyl_symbol_kernel_transport_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="weyl_symbol_kernel_transport_consequence_theorem.json",
    )
    args = parser.parse_args()

    theorem = load(args.transport_json)
    closed = bool(theorem.get("weylSymbolKernelTransportClosed"))

    transport_status = status(
        "Weyl symbol to coordinate-kernel transport consequence",
        closed,
        (
            "In the fixed Weyl convention, the coordinate kernel quadratic "
            "form is exactly the Weyl operator quadratic form "
            "<Op^W(sigma_omega)f,f> on the test domain."
        ),
        blocker=None if closed else "Close weyl_symbol_kernel_transport_theorem.py.",
    )

    data = {
        "theoremName": "Weyl symbol coordinate-kernel transport consequence theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "weylSymbolKernelTransportConsequenceStatus": transport_status,
        },
        "weylSymbolKernelTransportClosed": transport_status["closed"],
        "coordinateKernelQuadraticFormTransportClosed": transport_status["closed"],
        "proof": [
            "Import the full Weyl symbol coordinate-kernel transport theorem.",
            "Export only the terminal transport consequence needed by the external foundation.",
        ],
        "remainingAnalyticGap": None if transport_status["closed"] else transport_status["blocker"],
        "nextProofTarget": None if transport_status["closed"] else transport_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Weyl symbol coordinate-kernel transport consequence theorem")
    print(f"  transport consequence closed: {data['weylSymbolKernelTransportClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
