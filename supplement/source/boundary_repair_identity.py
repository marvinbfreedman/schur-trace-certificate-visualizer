#!/usr/bin/env python3
r"""Boundary-repair identity ledger.

The quotient-to-original lift left two nominal options:

  (A) every primitive original test F satisfies R_global F=0;
  (B) the original Weyl quadratic form is Q_Phi(F)+||S R_global F||^2.

This ledger proves that (A) is not viable in the natural compact-smooth
primitive class and explains why (B) needs a canonical repair operator.  The
abstract quotient theorem supplies an S, but that S is not unique: if

    Q = P - R^* D R,      P>=0, D>=0,

then for every T>=0 on the trace range,

    Q = (P+R^*TR) - R^*(D+T)R.

Thus a boundary-repair identity cannot use an arbitrary quotient S.  It must
identify the canonical trace-side boundary form coming from the Weyl/Volterra
Green/Lagrange identity and prove that it equals an admissible quotient repair.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quotient-json",
        default="weyl_volterra_quotient_schur_theorem.json",
    )
    parser.add_argument("--json-out", default="boundary_repair_identity.json")
    args = parser.parse_args()

    quotient = load(args.quotient_json)
    quotient_closed = bool(
        quotient.get("globalWeylVolterraSchurStatus", {}).get("closed")
    )

    primitive_vanishing = status(
        "primitive-vanishing route is false",
        True,
        (
            "For any active a, Lambda_a is a nonzero "
            "order-eight jet functional.  By the elementary jet-extension "
            "lemma, choose a compactly supported smooth primitive F with "
            "prescribed jet j_a^8 F and Lambda_a(F) != 0.  Then f=F' is an "
            "original compact smooth test with zero total integral, while its "
            "primitive F is not in ker R_global."
        ),
    )
    alternative_a_resolution = status(
        "alternative A cannot be proved",
        True,
        (
            "Alternative (A) asserts primitive original tests satisfy "
            "R_global F=0.  The jet-extension proof above gives a compact "
            "smooth primitive F with Lambda_a(F) != 0 for some active a.  "
            "Thus (A) is false, not merely unproved."
        ),
    )
    repair_nonunique = status(
        "abstract quotient repair is nonunique",
        True,
        (
            "If Q=P-R^*DR with P>=0 and D>=0, then for every positive trace "
            "operator T, Q=(P+R^*TR)-R^*(D+T)R.  Therefore the abstract S from "
            "the quotient theorem is not canonical and cannot be identified "
            "with the original Weyl boundary contribution without a separate "
            "Green/Lagrange boundary calculation."
        ),
    )
    alternative_b_as_stated = status(
        "alternative B as stated is not well-defined",
        True,
        (
            "Alternative (B) asks for the repair with 'the same S from the "
            "quotient Schur theorem.'  But the quotient theorem supplies a "
            "nonunique repair family: D can be replaced by D+T for any T>=0 "
            "if the positive part P is also replaced by P+R^*TR.  Therefore "
            "there is no distinguished 'same S' until the Weyl/Volterra "
            "Green identity constructs a canonical boundary operator D_bdy."
        ),
    )
    canonical_boundary_form = status(
        "canonical Weyl/Volterra boundary form constructed",
        False,
        (
            "The current notes contain the Lagrange concomitant identity "
            "D_s B_P[h,f]=h P f-f P^*h and source-row constructions, but they "
            "do not yet assemble the original primitive Weyl form into a "
            "positive interior term plus a canonical trace quadratic form "
            "B_bdy(RF)."
        ),
        blocker=(
            "Derive the exact primitive Weyl Green identity and isolate the "
            "boundary quadratic form B_bdy on the completed trace range."
        ),
    )
    canonical_equals_quotient = status(
        "canonical boundary form equals an admissible quotient repair",
        False,
        (
            "Even after B_bdy is constructed, one must prove "
            "B_bdy(RF)=||S R_global F||^2 for an S that is compatible with "
            "the Douglas/Moore-Penrose repair in the quotient Schur theorem, "
            "or at least B_bdy >= S^*S on the trace range."
        ),
        blocker=(
            "Compare the canonical boundary operator D_bdy with the quotient "
            "repair operator D_q=(Gamma^*Gamma-C)_+ on the transported trace "
            "range; prove equality or domination."
        ),
    )
    boundary_repair_identity = status(
        "boundary-repair identity",
        False,
        (
            "The identity original Weyl form = Q_Phi(F)+||S R_global F||^2 "
            "is not closed.  The primitive-vanishing alternative is false, "
            "and the abstract quotient repair is nonunique.  The remaining "
            "work is to construct the canonical boundary form and identify it "
            "with an admissible quotient repair."
        ),
        blocker=(
            "Construct D_bdy from the exact Weyl/Volterra Green identity and "
            "prove D_bdy = D_q, or D_bdy >= D_q, on the completed trace range."
        ),
    )

    data = {
        "theoremName": "boundary-repair identity",
        "quotientJson": args.quotient_json,
        "normalizedSchurCertificateClosed": quotient_closed,
        "primitiveVanishingFalse": True,
        "primitiveVanishingRouteStatus": primitive_vanishing,
        "alternativeAResolutionStatus": alternative_a_resolution,
        "abstractRepairNonuniquenessStatus": repair_nonunique,
        "alternativeBAsStatedStatus": alternative_b_as_stated,
        "canonicalBoundaryFormStatus": canonical_boundary_form,
        "canonicalEqualsQuotientRepairStatus": canonical_equals_quotient,
        "boundaryRepairIdentityStatus": boundary_repair_identity,
        "boundaryRepairIdentityClosed": False,
        "proofOfPrimitiveVanishingFailure": [
            (
                "For active a, the defect row e(a) is nonzero, hence "
                "Lambda_a is a nonzero continuous linear functional on the "
                "finite jet space J_a^8."
            ),
            (
                "Choose a jet v with <e(a),v> != 0.  By the compactly "
                "supported jet-extension lemma, there is F in C_c^infty with "
                "j_a^8F=v."
            ),
            (
                "Then f=F' is an admissible compact smooth original test "
                "with primitive F, but (R_global F)(a)=Lambda_a(F) != 0.  "
                "Therefore primitive original tests are not contained in "
                "ker R_global."
            ),
        ],
        "nonuniquenessProof": (
            "Given any valid repair Q=P-R^*DR and any bounded T>=0 on the "
            "trace range, set P_T=P+R^*TR and D_T=D+T.  Then "
            "P_T>=0, D_T>=0, and Q=P_T-R^*D_T R.  Hence the abstract quotient "
            "repair is a family, not a canonical boundary term."
        ),
        "requestedAlternativesConclusion": (
            "Neither requested alternative can be proved as stated: (A) is "
            "false, and (B) is not well-defined until a canonical boundary "
            "repair D_bdy is derived.  The corrected provable target is "
            "D_bdy = D_q, or the weaker sufficient comparison D_bdy >= D_q, "
            "on the completed trace range."
        ),
        "nextProofTarget": boundary_repair_identity["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Boundary-repair identity")
    print(f"  primitive-vanishing disproof closed: {primitive_vanishing['closed']}")
    print(f"  alternative A false: {alternative_a_resolution['closed']}")
    print(f"  abstract repair nonuniqueness proved: {repair_nonunique['closed']}")
    print(f"  alternative B as stated not well-defined: {alternative_b_as_stated['closed']}")
    print(f"  canonical boundary form constructed: {canonical_boundary_form['closed']}")
    print(f"  boundary repair identity closed: {boundary_repair_identity['closed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
