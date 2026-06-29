#!/usr/bin/env python3
r"""Finite PSD cone closure theorem.

This theorem records the elementary closure fact used in the augmented
closed-cone limit:

If G_n are N x N Hermitian positive semidefinite matrices and G_n converges
entrywise to G, then G is Hermitian positive semidefinite.

Proof: for every vector c in C^N,

    c^* G c = lim_n c^* G_n c >= 0.

This is a finite-dimensional symbolic theorem.  It has no numerical input.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", default="finite_psd_cone_closure_theorem.json")
    args = parser.parse_args()

    hermitian_status = status(
        "Hermitian entrywise limit",
        True,
        "Entrywise limits of Hermitian matrices remain Hermitian.",
    )
    quadratic_status = status(
        "quadratic-form limit",
        True,
        (
            "For fixed c in C^N, the scalar c^*G_n c is a finite linear "
            "combination of entries of G_n, so it converges to c^*G c."
        ),
    )
    closure_status = status(
        "finite PSD cone closed",
        True,
        (
            "Since c^*G_n c>=0 for every n and every c, the limit satisfies "
            "c^*G c>=0 for every c.  Hence G is positive semidefinite."
        ),
    )

    data = {
        "theoremName": "finite PSD cone closure theorem",
        "proofClass": "symbolic identity",
        "statuses": {
            "hermitianLimitStatus": hermitian_status,
            "quadraticFormLimitStatus": quadratic_status,
            "finitePsdConeClosureStatus": closure_status,
        },
        "finitePsdConeClosureClosed": True,
        "entrywiseLimitPreservesPsd": True,
        "proof": [
            "Fix a finite dimension N and a vector c in C^N.",
            "Entrywise convergence G_n->G implies c^*G_n c -> c^*G c.",
            "Each c^*G_n c is nonnegative because G_n is PSD.",
            "Therefore c^*G c>=0 for all c, which is exactly G>=0.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Finite PSD cone closure theorem")
    print(f"  closure closed: {closure_status['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
