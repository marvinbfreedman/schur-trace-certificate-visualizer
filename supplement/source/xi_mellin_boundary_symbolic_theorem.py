#!/usr/bin/env python3
r"""Symbolic Mellin boundary theorem for the augmented Xi trace.

This theorem object extracts the exact analytic identities from the older
numeric Mellin split and boundary-concomitant scripts.

For one finite Xi atom

    f(v)=a exp(beta v-c exp(v)),      v>=0,

and Fourier scale tau, define alpha=beta+i tau z.  Then

    X(z)=1/2 int_0^infty f(v) exp(i tau z v) dv
        =1/2 a c^{-alpha} Gamma(alpha,c).

Splitting at a Volterra base point s gives the exact decomposition

    X(z)=B(s,z)+T(s,z),

where

    B(s,z)=1/2 a c^{-alpha} int_c^{c exp(s)} y^{alpha-1} exp(-y) dy

and the tail is the diagonal Volterra ratio atom after the substitution
v=s+u.  The prefix satisfies

    d_s B(s,z)=1/2 a exp(beta s-c exp(s)) exp(i tau z s),
    B(0,z)=0.

Thus the missing boundary coordinate is the primitive trace

    Mu_z(f)=int_0^L B(s,z)f(s)ds

or, with F(s)=int_s^L f(t)dt,

    Mu_z(f)=int_0^L B'(s,z)F(s)ds,

because B(0,z)=0 and F(L)=0.
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
    parser.add_argument("--json-out", default="xi_mellin_boundary_symbolic_theorem.json")
    args = parser.parse_args()

    mellin_total = status(
        "Xi atom Mellin integral",
        True,
        (
            "With y=c exp(v), dv=dy/y and exp(beta v)exp(i tau z v) "
            "=(y/c)^(beta+i tau z).  This gives "
            "X(z)=1/2 a c^{-alpha} Gamma(alpha,c)."
        ),
    )
    split = status(
        "moving boundary split",
        True,
        (
            "Splitting the y-integral at c exp(s) gives the incomplete-gamma "
            "prefix B(s,z) plus the tail integral from c exp(s) to infinity."
        ),
    )
    volterra_tail = status(
        "diagonal Volterra tail transport",
        True,
        (
            "In the tail use v=s+u.  The factor exp(beta(s+u)-c exp(s+u)) "
            "separates into the base value at s times "
            "exp(beta u-c exp(s)(exp(u)-1)), exactly the Volterra ratio atom."
        ),
    )
    boundary = status(
        "Mellin boundary concomitant",
        True,
        (
            "Differentiating the prefix endpoint c exp(s) gives "
            "B'(s,z)=1/2 a exp(beta s-c exp(s))exp(i tau z s), and B(0,z)=0."
        ),
    )
    primitive = status(
        "primitive Mu trace",
        True,
        (
            "For F(s)=int_s^L f(t)dt, integration by parts gives "
            "int_0^L B(s,z)f(s)ds=int_0^L B'(s,z)F(s)ds because "
            "B(0,z)=0 and F(L)=0."
        ),
    )

    data = {
        "theoremName": "symbolic Mellin boundary theorem",
        "proofClass": "symbolic identity",
        "symbols": {
            "atom": "f(v)=a exp(beta v-c exp(v)), v>=0",
            "alpha": "alpha=beta+i tau z",
            "total": "X(z)=1/2 a c^{-alpha} Gamma(alpha,c)",
            "prefix": (
                "B(s,z)=1/2 a c^{-alpha} int_c^{c exp(s)} "
                "y^{alpha-1} exp(-y) dy"
            ),
            "tail": "T(s,z)=X(z)-B(s,z)",
            "mu": "Mu_z(f)=int_0^L B(s,z)f(s)ds",
        },
        "statuses": {
            "mellinTotalStatus": mellin_total,
            "movingBoundarySplitStatus": split,
            "diagonalVolterraTailStatus": volterra_tail,
            "boundaryConcomitantStatus": boundary,
            "primitiveMuTraceStatus": primitive,
        },
        "exactMellinSplitClosed": True,
        "mellinBoundaryConcomitantClosed": True,
        "primitiveMuTraceClosed": True,
        "finiteAtomTransportClosed": True,
        "proof": [
            "Substitute y=c exp(v) in the half-line Mellin atom integral.",
            "Split the resulting incomplete-gamma integral at y=c exp(s).",
            "Use v=s+u in the tail to identify the normalized Volterra ratio atom.",
            "Differentiate the moving upper endpoint to obtain B'(s,z).",
            "Integrate by parts against F(s)=int_s^L f(t)dt to get the Mu trace.",
        ],
        "scope": (
            "This theorem proves only the exact symbolic atom identities.  "
            "Trace closure, repair positivity, and positive-cone limits are "
            "handled by separate theorem objects."
        ),
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Symbolic Mellin boundary theorem")
    print(f"  exact Mellin split: {data['exactMellinSplitClosed']}")
    print(f"  boundary concomitant: {data['mellinBoundaryConcomitantClosed']}")
    print(f"  primitive Mu trace: {data['primitiveMuTraceClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
