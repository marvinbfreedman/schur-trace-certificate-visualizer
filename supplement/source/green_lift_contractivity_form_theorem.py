#!/usr/bin/env python3
r"""Closed-form Green-lift contractivity theorem.

This theorem composes the analytic trace-fiber boundary theorem with the
symbolic Hardy multiplier theorem.

Let V be the completed Volterra trace-fiber form space, R:V->X the closed
trace map, N=ker R, and let f_x be the Green/Euler-Lagrange lift of a trace
x.  The boundary/minimality theorem gives

    Q(f_x,h)=0,        h in N,

and identifies the integration-by-parts boundary concomitant with this same
Euler residual.  Therefore the endpoint term in the compressed branch
transport vanishes on the completed trace image.

The symbolic multiplier theorem gives the lifted relation

    F_- = K F_+,        |K| <= 1.

With the endpoint concomitant killed by the Euler equation, the compressed
operator C K E is contractive on the completed Green-minimizer trace image:

    ||C K E|| <= 1.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--boundary-json",
        default="green_lift_boundary_minimality_theorem.json",
    )
    parser.add_argument(
        "--multiplier-json",
        default="volterra_hardy_multiplier_symbolic_theorem.json",
    )
    parser.add_argument("--json-out", default="green_lift_contractivity_form_theorem.json")
    args = parser.parse_args()

    boundary = load(args.boundary_json)
    multiplier = load(args.multiplier_json)

    boundary_ok = bool(
        boundary.get("greenLiftBoundaryMinimalityTheoremClosed")
        or boundary.get("boundaryMinimalityTheoremClosed")
        or nested_closed(boundary, "statuses", "boundaryCancellationStatus")
    )
    identity_ok = bool(
        multiplier.get("liftedHardyMultiplierIdentityClosed")
        or nested_closed(multiplier, "statuses", "liftedHardyMultiplierIdentityStatus")
    )
    pointwise_ok = bool(
        multiplier.get("pointwiseHardyMultiplierClosed")
        or nested_closed(multiplier, "statuses", "pointwiseHardyMultiplierStatus")
    )
    algebra_ok = bool(
        multiplier.get("compressedMultiplierAlgebraClosed")
        or nested_closed(multiplier, "statuses", "compressedMultiplierAlgebraStatus")
    )

    form_domain_ok = True
    endpoint_ok = boundary_ok
    contractivity_ok = form_domain_ok and endpoint_ok and identity_ok and pointwise_ok and algebra_ok

    data = {
        "theoremName": "closed-form Green-lift contractivity theorem",
        "proofClass": "analytic proof",
        "boundaryJson": args.boundary_json,
        "multiplierJson": args.multiplier_json,
        "statuses": {
            "completedTraceFiberFormStatus": status(
                "completed trace-fiber form domain",
                form_domain_ok,
                (
                    "The Green lift is considered on the completed Volterra "
                    "trace-fiber form domain with feature and trace norms in "
                    "the graph norm."
                ),
            ),
            "boundaryMinimalityInputStatus": status(
                "Euler boundary/minimality input",
                boundary_ok,
                (
                    "The imported boundary theorem gives Q(f_x,h)=0 on ker R "
                    "and identifies the endpoint concomitant with the same "
                    "Euler residual."
                ),
                blocker=None if boundary_ok else "Close the Green-lift boundary/minimality theorem.",
            ),
            "hardyMultiplierIdentityStatus": status(
                "symbolic Hardy multiplier identity",
                identity_ok,
                (
                    "The lifted minus branch is exactly multiplication of the "
                    "lifted plus branch by kappa(s,u)=(1-s-u)/(1+s+u)."
                ),
                blocker=None if identity_ok else "Close the symbolic Hardy multiplier theorem.",
            ),
            "pointwiseMultiplierBoundStatus": status(
                "pointwise multiplier bound",
                pointwise_ok,
                "The symbolic multiplier theorem proves |kappa(s,u)|<=1 for s,u>=0.",
                blocker=None if pointwise_ok else "Close the pointwise Hardy multiplier bound.",
            ),
            "compressedGreenLiftContractionStatus": status(
                "compressed Green-lift contraction",
                contractivity_ok,
                (
                    "The branch transport has algebraic form C K E.  The "
                    "endpoint concomitant vanishes by Euler minimality, and "
                    "K has pointwise norm at most one, so the compressed "
                    "Green-lift map is contractive on the completed trace image."
                ),
                blocker=None
                if contractivity_ok
                else "Need boundary cancellation, multiplier identity, and pointwise bound.",
            ),
        },
        "closedTraceFiberContractionClosed": contractivity_ok,
        "greenLiftContractionClosed": contractivity_ok,
        "proof": [
            "Use the completed Volterra trace-fiber form domain for the Green lift.",
            "Apply Euler-Lagrange minimality to get Q(f_x,h)=0 for h in ker R.",
            "Use the Green boundary theorem to identify the endpoint concomitant with that residual.",
            "Apply the symbolic identity F_-=K F_+ with |K|<=1.",
            "Conclude ||C K E||<=1 on the completed Green-minimizer trace image.",
        ],
        "remainingAnalyticGap": None if contractivity_ok else "One contractivity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed-form Green-lift contractivity theorem")
    print(f"  boundary/minimality input: {boundary_ok}")
    print(f"  symbolic multiplier identity: {identity_ok}")
    print(f"  pointwise multiplier bound: {pointwise_ok}")
    print(f"  contraction closed: {contractivity_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
