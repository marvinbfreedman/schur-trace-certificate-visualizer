#!/usr/bin/env python3
r"""Slim consequence for the trace Lagrange adjoint identity.

The canonical boundary repair comparison only needs to know that the moving
trace Lagrange identity is exact.  Older ledgers consumed this as a numerical
``maxIdentityRelativeDefect`` field from ``trace_lagrange_adjoint_control``.
This wrapper exports the same compatibility key from the symbolic identity
theorem, with defect zero.
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
        "--identity-json",
        default="trace_lagrange_adjoint_identity_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_lagrange_adjoint_identity_consequence_theorem.json",
    )
    args = parser.parse_args()

    identity = load(args.identity_json)
    identity_status = identity.get("lagrangeIdentityStatus", {})
    ok = bool(identity.get("theoremClosed")) and bool(identity_status.get("closed"))

    data = {
        "theoremName": "trace Lagrange adjoint identity consequence theorem",
        "proofClass": "symbolic identity",
        "identityJson": args.identity_json,
        "statuses": {
            "lagrangeIdentityStatus": status(
                "exact moving-trace Lagrange identity",
                ok,
                (
                    "The concomitant identity D_s B_P[h,f]=hPf-fP^*h is an "
                    "algebraic product-rule telescoping identity for arbitrary "
                    "smooth trace coefficients a_k(s)."
                ),
            )
        },
        "lagrangeIdentityStatus": status(
            "exact moving-trace Lagrange identity",
            ok,
            (
                "The symbolic trace Lagrange adjoint identity theorem proves "
                "D_s B_P[h,f]=hPf-fP^*h exactly."
            ),
        ),
        "theoremClosed": ok,
        "lagrangeAdjointIdentityClosed": ok,
        "maxIdentityDefect": 0.0 if ok else 1.0,
        "maxIdentityRelativeDefect": 0.0 if ok else 1.0,
        "proof": [
            "Import the symbolic trace Lagrange adjoint identity theorem.",
            "Expose the legacy zero-defect scalar consumed by older comparison ledgers.",
            "No quadrature or sampled trace rows enter this consequence.",
        ],
        "remainingAnalyticGap": None if ok else "Symbolic trace Lagrange identity theorem is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace Lagrange adjoint identity consequence theorem")
    print(f"  Lagrange identity: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
