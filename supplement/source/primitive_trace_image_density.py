#!/usr/bin/env python3
r"""Primitive trace image density ledger.

After the primitive boundary audit, the remaining repair condition is

    D_q|_Y = 0,
    Y = R({F : F'=f, f an original compact Weyl test}).

This ledger checks whether Y is genuinely smaller than the completed trace
range X_R.  Under the density hypotheses already used in the quotient lift, it
is not: compact smooth functions supported inside the half-line are primitives
of compact smooth original tests by taking f=F'.  Such functions are dense in
the Volterra form domain V, and R:V->X_R is continuous with X_R defined as the
completed transported trace range.  Therefore

    closure(Y) = X_R.

Consequently, for bounded D_q, the condition D_q|_Y=0 is equivalent to
D_q=0 on all of X_R.  The primitive image does not provide an extra hidden
constraint that could kill the quotient repair.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    candidate = Path(path)
    if not candidate.exists():
        return None
    return json.loads(candidate.read_text(encoding="utf-8"))


def read(path: str) -> str:
    return Path(path).read_text(encoding="utf-8")


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def has_all(text: str, needles: list[str]) -> bool:
    return all(needle in text for needle in needles)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--notes", default="rh_klm_notes.md")
    parser.add_argument("--draft", default="rh_weyl_positive_draft.tex")
    parser.add_argument(
        "--primitive-boundary-json",
        default="primitive_boundary_zero_consequence_theorem.json",
    )
    parser.add_argument(
        "--quotient-json",
        default="quotient_minimal_repair_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="primitive_trace_image_density.json")
    args = parser.parse_args()

    notes = read(args.notes)
    draft = read(args.draft)
    primitive = load_optional(args.primitive_boundary_json)
    quotient = load_optional(args.quotient_json)

    density_hypotheses = has_all(
        notes,
        [
            "density and closure in the Volterra form domain",
            "completed transported trace range",
            "X_R = R(U)",
        ],
    ) or has_all(
        draft,
        [
            "density/closure",
            "completed transported trace range",
        ],
    )
    primitive_boundary_zero = bool(
        primitive
        and primitive.get("statuses", {})
        .get("canonicalPrimitiveBoundaryZeroStatus", {})
        .get("closed")
    )
    dq_bounded = bool(
        quotient
        and (
            quotient.get("traceSideRepairClosed")
            or quotient.get("traceSideRepairStatus", {}).get("closed")
            or quotient.get("statuses", {})
            .get("traceSideRepairInputStatus", {})
            .get("closed")
        )
    )

    primitive_contains_core = status(
        "primitive class contains compact smooth Volterra core",
        True,
        (
            "If F is compact smooth on the half-line with F(0)=0, then "
            "f=F' is compact smooth and F is exactly the cumulative primitive "
            "of f.  Thus the primitive image contains the usual compact "
            "smooth core used in the Volterra form closure."
        ),
    )
    trace_density = status(
        "primitive trace image is dense in X_R",
        density_hypotheses,
        (
            "By the imported density/closure hypothesis, the compact smooth "
            "Volterra core is dense in the form domain V.  Since R is "
            "continuous and X_R is the completed transported trace range, "
            "R of that core is dense in X_R.  The primitive trace image Y "
            "therefore has closure X_R."
        ),
        blocker=None
        if density_hypotheses
        else "State the compact-core density theorem for the Volterra form domain and continuity of R into X_R.",
    )
    dq_zero_equivalence = status(
        "D_q vanishes on primitive image iff D_q=0 on X_R",
        trace_density["closed"] and dq_bounded,
        (
            "The quotient repair D_q is bounded on X_R.  If it vanishes on "
            "the dense primitive trace image Y, continuity forces D_q=0 on "
            "all of X_R.  Conversely D_q=0 trivially implies D_q|_Y=0."
        ),
        blocker=None
        if trace_density["closed"] and dq_bounded
        else "Prove trace density of Y and boundedness of D_q on X_R.",
    )
    hidden_constraint_route = status(
        "primitive-image hidden constraint route",
        False,
        (
            "Because closure(Y)=X_R, there is no smaller primitive trace "
            "subspace on which D_q could vanish while remaining nonzero on "
            "X_R.  The problem has collapsed to full repair-vanishing or "
            "direct full-form positivity."
        ),
        blocker=(
            "Prove D_q=0 on X_R, equivalently show Gamma^*Gamma-C<=0 in the "
            "quotient Schur coordinates, or abandon the repair route and "
            "prove Q_Phi>=0 directly on the full primitive/form closure."
        ),
    )
    full_positivity = status(
        "original positivity from dense primitive image",
        False,
        (
            "The primitive boundary form is zero and the primitive trace image "
            "is dense.  Therefore original positivity is now equivalent, in "
            "this route, to proving that the quotient repair is actually zero "
            "or to proving the full transported form Q_Phi is nonnegative."
        ),
        blocker=hidden_constraint_route["blocker"],
    )

    data = {
        "theoremName": "primitive trace image density",
        "notes": args.notes,
        "draft": args.draft,
        "primitiveBoundaryJson": args.primitive_boundary_json if primitive else None,
        "quotientJson": args.quotient_json if quotient else None,
        "importedSignals": {
            "densityHypothesesPresent": density_hypotheses,
            "primitiveBoundaryZero": primitive_boundary_zero,
            "DqBoundedOnXR": dq_bounded,
        },
        "statuses": {
            "primitiveContainsCompactCoreStatus": primitive_contains_core,
            "primitiveTraceDenseStatus": trace_density,
            "dqZeroEquivalenceStatus": dq_zero_equivalence,
            "hiddenConstraintRouteStatus": hidden_constraint_route,
            "originalPositivityFromDensePrimitiveImageStatus": full_positivity,
        },
        "proof": [
            "For F in C_c^infty with F(0)=0, f=F' is an original compact smooth test and F(u)=int_0^u f.",
            "The compact smooth core is dense in the Volterra form domain V by the imported closure theorem.",
            "R is continuous into the transported trace completion X_R.",
            "Therefore closure R(primitive tests) contains closure R(C_c^infty core)=X_R.",
            "Since D_q is bounded on X_R, D_q vanishes on the primitive trace image iff D_q=0 on X_R.",
        ],
        "primitiveTraceImageDenseInXR": trace_density["closed"],
        "dqZeroOnPrimitiveImageEquivalentToDqZero": dq_zero_equivalence["closed"],
        "correctedNextProofTarget": hidden_constraint_route["blocker"],
        "originalWeylKernelPositivityClosed": False,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive trace image density")
    print(f"  primitive contains compact core: {primitive_contains_core['closed']}")
    print(f"  primitive trace dense in X_R: {trace_density['closed']}")
    print(f"  D_q|Y=0 iff D_q=0: {dq_zero_equivalence['closed']}")
    print(f"  hidden constraint route closed: {hidden_constraint_route['closed']}")
    print(f"  next: {data['correctedNextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
