#!/usr/bin/env python3
r"""Slim consequence for the quotient/original transport identity.

The publication quotient-to-original Weyl lift theorem only needs the closed
transport identity flag: parity reduction, primitive mixed-derivative
transport, Volterra/log normalization, and form-core closure have already been
assembled into a single exact identity theorem.
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
        "--transport-json",
        default="quotient_original_transport_identity_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="quotient_original_transport_identity_consequence_theorem.json",
    )
    args = parser.parse_args()

    transport = load(args.transport_json)
    ok = bool(transport.get("quotientOriginalTransportIdentityClosed"))

    data = {
        "theoremName": "quotient/original transport identity consequence theorem",
        "proofClass": "symbolic identity",
        "transportSource": "quotient original transport identity theorem",
        "statuses": {
            "transportIdentityStatus": status(
                "quotient/original transport identity",
                ok,
                (
                    "The imported transport identity theorem closes parity "
                    "reduction, primitive mixed-derivative transport, "
                    "Volterra/log normalization, and form-core closure."
                ),
            )
        },
        "quotientOriginalTransportIdentityClosed": ok,
        "proof": [
            "Import the closed quotient/original transport identity theorem.",
            "Expose only the transport-identity flag consumed by the publication lift theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Transport identity theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Quotient/original transport identity consequence theorem")
    print(f"  transport identity: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
