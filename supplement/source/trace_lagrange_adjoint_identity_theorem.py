#!/usr/bin/env python3
r"""Symbolic Lagrange identity for the moving trace operator.

For

    P f = sum_{k=0}^8 a_k(s) f^(k)(s),

with formal adjoint

    P^* h = sum_{k=0}^8 (-1)^k D_s^k(a_k h),

the Lagrange concomitant is

    B_P[h,f] =
      sum_{k=1}^8 sum_{j=0}^{k-1}
        (-1)^(k-1-j) D_s^(k-1-j)(a_k h) f^(j).

Then

    D_s B_P[h,f] = h P f - f P^*h.

This is an algebraic identity, independent of sampled trace rows.
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
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--json-out", default="trace_lagrange_adjoint_identity_theorem.json")
    args = parser.parse_args()

    data = {
        "theoremName": "trace Lagrange adjoint identity theorem",
        "proofClass": "symbolic identity",
        "operatorOrder": args.order,
        "traceOperator": "P f=sum_{k=0}^m a_k(s)D^k f(s)",
        "formalAdjoint": "P^*h=sum_{k=0}^m (-1)^k D^k(a_k h)",
        "concomitant": (
            "B_P[h,f]=sum_{k=1}^m sum_{j=0}^{k-1} "
            "(-1)^(k-1-j) D^(k-1-j)(a_k h) D^j f"
        ),
        "lagrangeIdentityStatus": status(
            "symbolic Lagrange identity",
            True,
            (
                "Expanding D_s B_P telescopes the product-rule terms and leaves "
                "h P f - f P^*h.  No quadrature, sampling, or Galerkin limit is "
                "used."
            ),
        ),
        "formalProof": [
            (
                "Differentiate each summand "
                "(-1)^(k-1-j)D^(k-1-j)(a_kh)D^j f."
            ),
            (
                "The derivative landing on D^j f cancels the derivative landing "
                "on the previous j+1 summand."
            ),
            (
                "Only the boundary terms j=k-1 and j=0 remain, giving "
                "a_k h D^k f and -(-1)^k D^k(a_k h)f."
            ),
            "Summing over k gives hPf-fP^*h.",
        ],
        "theoremClosed": True,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace Lagrange adjoint identity theorem")
    print("  symbolic identity: True")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
