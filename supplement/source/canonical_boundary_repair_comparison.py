#!/usr/bin/env python3
r"""Canonical boundary repair comparison ledger.

The quotient Schur theorem gives a noncanonical factorization

    Q_Phi(f) = ||G_q f||^2 - <D_q Rf, Rf>,

where

    D_q = (Gamma^* Gamma - C)_+

is the minimal Moore--Penrose/Douglas repair on the transported trace range
when the positive part is fixed to be the quotient Schur positive form.

The original Weyl lift cannot use an arbitrary repair.  The repair must be the
canonical boundary form produced by the exact primitive Weyl/Volterra Green
identity.  If

    beta_bdy(f,g) := Q_original(f,g) - Q_Phi(f,g)

descends through R, then it defines a unique self-adjoint trace operator
D_bdy by

    beta_bdy(f,g) = <D_bdy Rf, Rg>_{X_R}.

This script records the abstract theorem:

    D_bdy >= D_q  =>  Q_original(f) >= 0.

It does not mark the analytic comparison closed until an exact Green-boundary
calculation proves the descent of beta_bdy and the operator inequality
D_bdy-D_q >= 0 on the completed trace range.
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


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def nested_closed(data: dict | None, key: str) -> bool:
    if not data:
        return False
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quotient-json",
        default="quotient_minimal_repair_consequence_theorem.json",
    )
    parser.add_argument(
        "--boundary-repair-json",
        default="boundary_repair_diagnostic_consequence_theorem.json",
    )
    parser.add_argument(
        "--lagrange-json",
        default="trace_lagrange_adjoint_identity_consequence_theorem.json",
    )
    parser.add_argument(
        "--trace-green-json",
        default="continuum_trace_to_source_green_diagnostic_consequence_theorem.json",
    )
    parser.add_argument(
        "--primitive-boundary-json",
        default="primitive_boundary_zero_consequence_theorem.json",
    )
    parser.add_argument(
        "--primitive-density-json",
        default="primitive_trace_density_consequence_theorem.json",
    )
    parser.add_argument(
        "--primitive-endpoint-json",
        default="primitive_endpoint_compatibility_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="canonical_boundary_repair_comparison.json",
    )
    args = parser.parse_args()

    quotient = load_optional(args.quotient_json)
    repair = load_optional(args.boundary_repair_json)
    lagrange = load_optional(args.lagrange_json)
    trace_green = load_optional(args.trace_green_json)
    primitive_boundary = load_optional(args.primitive_boundary_json)
    primitive_density = load_optional(args.primitive_density_json)
    primitive_endpoint = load_optional(args.primitive_endpoint_json)

    quotient_closed = bool(
        (quotient or {}).get("globalWeylVolterraSchurClosed")
        or nested_closed(quotient, "globalWeylVolterraSchurStatus")
    )
    quotient_factorization_closed = bool(
        (quotient or {}).get("quotientFactorizationClosed")
        or nested_closed(quotient, "quotientFactorizationStatus")
    )
    trace_repair_closed = bool(
        (quotient or {}).get("traceSideRepairClosed")
        or nested_closed(quotient, "traceSideRepairStatus")
    )
    minimal_repair_closed = bool(
        (quotient or {}).get("minimalQuotientRepairClosed")
        or (quotient or {}).get("quotientMinimalRepairConsequenceClosed")
    )
    primitive_route_false = bool(repair and repair.get("primitiveVanishingFalse"))
    abstract_repair_nonunique = nested_closed(repair, "abstractRepairNonuniquenessStatus")
    lagrange_identity_closed = (
        lagrange is not None
        and float(lagrange.get("maxIdentityRelativeDefect", 1.0)) < 1e-50
    )

    minimal_repair = status(
        "minimal quotient repair D_q defined",
        minimal_repair_closed
        or (quotient_closed and quotient_factorization_closed and trace_repair_closed),
        (
            "The quotient theorem fixes the minimal trace repair for its "
            "chosen positive part: D_q=(Gamma^*Gamma-C)_+ on the completed "
            "transported trace range X_R.  This is invariant only after the "
            "positive quotient Schur form has been fixed."
        ),
        blocker=None
        if quotient_closed
        else "Close the normalized full-Phi quotient Schur certificate.",
    )
    canonical_needed = status(
        "canonical repair is necessary",
        primitive_route_false and abstract_repair_nonunique,
        (
            "Primitive original tests are not contained in ker R_global, and "
            "the abstract quotient repair is nonunique.  Therefore the "
            "original lift must use the Green-boundary repair D_bdy rather "
            "than an arbitrary S from the quotient factorization."
        ),
        blocker=None
        if primitive_route_false and abstract_repair_nonunique
        else "Retire the false primitive-vanishing route and the noncanonical S route.",
    )
    abstract_descent = status(
        "abstract boundary descent theorem",
        True,
        (
            "Let beta(f,g)=Q_original(f,g)-Q_Phi(f,g).  Beta descends to a "
            "unique bounded Hermitian form on X_R exactly when beta(n,f)=0 "
            "for every n in ker R and beta is continuous in the transported "
            "trace norm.  Then D_bdy is defined by "
            "<D_bdy Rf,Rg>_X_R=beta(f,g)."
        ),
    )
    comparison_sufficiency = status(
        "comparison implies original positivity",
        True,
        (
            "If Q_Phi=||G_q f||^2-<D_q Rf,Rf>, beta=<D_bdy Rf,Rf>, and "
            "D_bdy>=D_q, then Q_original=||G_q f||^2+"
            "<(D_bdy-D_q)Rf,Rf> >= 0.  Equality D_bdy=D_q is the sharp "
            "boundary-repair identity."
        ),
    )
    green_identity_input = status(
        "exact Green identity input available",
        lagrange_identity_closed,
        (
            "The moving-trace Lagrange identity "
            "D_s B_P[h,f]=hP f-fP^*h is verified to high precision and is "
            "the local identity needed to derive endpoint trace terms."
        ),
        blocker=None
        if lagrange_identity_closed
        else "Verify the moving-trace Lagrange identity and source concomitant.",
    )

    primitive_boundary_zero = bool(
        primitive_boundary
        and primitive_boundary.get("statuses", {})
        .get("canonicalPrimitiveBoundaryZeroStatus", {})
        .get("closed")
    )
    primitive_boundary_descends = bool(
        primitive_boundary
        and primitive_boundary.get("statuses", {})
        .get("zeroBoundaryDescendsStatus", {})
        .get("closed")
    )
    primitive_trace_dense = bool(
        primitive_density
        and primitive_density.get("primitiveTraceImageDenseInXR")
    )
    dq_zero_equivalent = bool(
        primitive_density
        and primitive_density.get("dqZeroOnPrimitiveImageEquivalentToDqZero")
    )
    primitive_endpoint_closed = bool(
        primitive_endpoint
        and primitive_endpoint.get("primitiveEndpointCompatibilityClosed")
    )
    primitive_endpoint_dq_zero = bool(
        primitive_endpoint
        and primitive_endpoint.get("statuses", {})
        .get("dqVanishesOnXRStatus", {})
        .get("closed")
    )

    # The primitive transport audit shows that the x,y Green endpoints do not
    # generate a hidden positive repair; D_bdy is zero for the already
    # identified Q_original -> Q_Phi transport.
    trace_green_summary = (trace_green or {}).get("statuses", {})
    exact_primitive_boundary = status(
        "canonical D_bdy constructed from primitive Weyl Green identity",
        primitive_boundary_zero,
        (
            (
                "The primitive transport identity has no lower or upper "
                "endpoint contribution: F(0)=0 kills the lower endpoint and "
                "kernel decay kills infinity.  Therefore "
                "beta_bdy=Q_original-Q_Phi=0 and D_bdy=0."
            )
            if primitive_boundary_zero
            else (
                "The current Green ledgers construct the active trace-to-source "
                "BVP and verify the local Lagrange identity.  They do not yet "
                "assemble Q_original-Q_Phi into an explicit completed-trace "
                "quadratic form beta_bdy(RF)."
            )
        ),
        blocker=None
        if primitive_boundary_zero
        else (
            "Derive the primitive Weyl/Volterra Green identity with all "
            "integration endpoints retained, then collect exactly the terms "
            "depending only on R_global F."
        ),
    )
    boundary_descends = status(
        "Green boundary form descends through R_global",
        primitive_boundary_descends,
        (
            (
                "Since beta_bdy is the zero form, it annihilates ker R_global "
                "and is bounded in the transported trace norm.  Thus it "
                "descends trivially with D_bdy=0."
            )
            if primitive_boundary_descends
            else (
                "After beta_bdy is written down, one must prove beta_bdy(n,f)=0 "
                "for every n in ker R_global and prove boundedness in the "
                "transported trace norm X_R.  This is the well-definedness of "
                "D_bdy on the completed trace range."
            )
        ),
        blocker=None
        if primitive_boundary_descends
        else (
            "Show the exact boundary terms are functions only of the moving "
            "trace Lambda_a(F), not of an arbitrary representative F."
        ),
    )
    domination = status(
        "canonical boundary dominates minimal quotient repair",
        primitive_endpoint_dq_zero,
        (
            (
                "The primitive endpoint compatibility theorem proves D_q=0 "
                "on X_R.  Since the primitive transport gives D_bdy=0, the "
                "comparison D_bdy>=D_q holds with equality."
            )
            if primitive_endpoint_dq_zero
            else (
                "The primitive transport gives D_bdy=0.  Therefore the "
                "comparison D_bdy>=D_q is equivalent to D_q=0 on the relevant "
                "primitive trace image.  This has not been proved."
            )
            if primitive_boundary_zero
            else (
                "The sufficient inequality is D_bdy-D_q>=0 on X_R.  No current "
                "certificate contains the matrix/kernel of D_bdy, so this "
                "operator comparison has not yet been tested or proved."
            )
        ),
        blocker=None
        if primitive_endpoint_dq_zero
        else (
            (
                (
                    "Prove D_q=0 on X_R, equivalently show "
                    "Gamma^*Gamma-C<=0 in the quotient Schur coordinates, or "
                    "prove Q_Phi>=0 directly on the full primitive/form closure."
                )
                if dq_zero_equivalent
                else (
                    "Characterize the primitive trace image Y=R({F: F'=original "
                    "test}) inside X_R and prove D_q|_Y=0, or prove the "
                    "equivalent direct positivity of Q_Phi on the primitive image."
                )
            )
            if primitive_boundary_zero
            else (
                "Represent D_bdy in the same transported trace coordinates as "
                "D_q=(Gamma^*Gamma-C)_+, then prove the positive semidefinite "
                "comparison D_bdy-D_q>=0."
            )
        ),
    )

    comparison_closed = (
        minimal_repair["closed"]
        and canonical_needed["closed"]
        and exact_primitive_boundary["closed"]
        and boundary_descends["closed"]
        and domination["closed"]
    )
    original_positivity = status(
        "original Weyl form positivity from boundary comparison",
        comparison_closed,
        (
            "When the canonical boundary form is constructed, descends to "
            "X_R, and dominates D_q, the original Weyl quadratic form is a "
            "sum of two nonnegative terms."
        ),
        blocker=None
        if comparison_closed
        else "Close the canonical D_bdy construction/descent and prove D_bdy>=D_q.",
    )

    data = {
        "theoremName": "canonical boundary repair comparison",
        "quotientJson": args.quotient_json if quotient else None,
        "boundaryRepairJson": args.boundary_repair_json if repair else None,
        "lagrangeJson": args.lagrange_json if lagrange else None,
        "traceGreenJson": args.trace_green_json if trace_green else None,
        "primitiveBoundaryJson": args.primitive_boundary_json if primitive_boundary else None,
        "primitiveDensityJson": args.primitive_density_json if primitive_density else None,
        "primitiveEndpointJson": args.primitive_endpoint_json if primitive_endpoint else None,
        "importedSignals": {
            "normalizedQuotientSchurClosed": quotient_closed,
            "quotientFactorizationClosed": quotient_factorization_closed,
            "traceSideRepairClosed": trace_repair_closed,
            "minimalQuotientRepairClosed": minimal_repair_closed,
            "primitiveVanishingRouteFalse": primitive_route_false,
            "abstractRepairNonunique": abstract_repair_nonunique,
            "lagrangeIdentityClosed": lagrange_identity_closed,
            "lagrangeMaxRelativeDefect": (lagrange or {}).get(
                "maxIdentityRelativeDefect"
            ),
            "traceGreenStatuses": trace_green_summary,
            "primitiveBoundaryZero": primitive_boundary_zero,
            "primitiveBoundaryDescends": primitive_boundary_descends,
            "primitiveTraceDenseInXR": primitive_trace_dense,
            "DqZeroOnPrimitiveImageEquivalentToDqZero": dq_zero_equivalent,
            "primitiveEndpointCompatibilityClosed": primitive_endpoint_closed,
            "primitiveEndpointDqZero": primitive_endpoint_dq_zero,
        },
        "statuses": {
            "minimalQuotientRepairStatus": minimal_repair,
            "canonicalRepairNeededStatus": canonical_needed,
            "abstractBoundaryDescentStatus": abstract_descent,
            "comparisonSufficiencyStatus": comparison_sufficiency,
            "greenIdentityInputStatus": green_identity_input,
            "canonicalBoundaryConstructionStatus": exact_primitive_boundary,
            "boundaryDescentThroughTraceStatus": boundary_descends,
            "boundaryDominatesQuotientRepairStatus": domination,
            "originalWeylPositivityFromComparisonStatus": original_positivity,
        },
        "operators": {
            "Dq": {
                "definition": "D_q=(Gamma^*Gamma-C)_+ on X_R=R(U) with ||Ru||_X_R=||u||_V.",
                "source": "closed quotient Schur theorem",
                "minimalForFixedPositivePart": True,
            },
            "Dbdy": {
                "definition": (
                    "If beta_bdy(f,g)=Q_original(f,g)-Q_Phi(f,g) descends "
                    "through R, define <D_bdy Rf,Rg>_X_R=beta_bdy(f,g)."
                ),
                "source": "exact primitive Weyl/Volterra Green boundary terms",
                "constructedHere": primitive_boundary_zero,
                "value": "0" if primitive_boundary_zero else "not constructed",
            },
        },
        "closedAbstractProof": [
            (
                "The descent criterion is exact: beta defines a trace-side "
                "operator iff beta vanishes whenever one argument lies in "
                "ker R and is bounded in the transported trace norm."
            ),
            (
                "With D_bdy defined by descent, "
                "Q_original=Q_Phi+<D_bdy Rf,Rf>."
            ),
            (
                "The quotient theorem gives "
                "Q_Phi=||G_q f||^2-<D_q Rf,Rf>."
            ),
            (
                "Therefore "
                "Q_original=||G_q f||^2+<((D_bdy-D_q)Rf),Rf>; "
                "D_bdy>=D_q proves original positivity."
            ),
        ],
        "canonicalBoundaryComparisonClosed": comparison_closed,
        "originalWeylKernelPositivityClosed": original_positivity["closed"],
        "nextProofTarget": (
            (
                (
                    "Prove D_q=0 on X_R, equivalently show "
                    "Gamma^*Gamma-C<=0 in the quotient Schur coordinates, or "
                    "prove Q_Phi>=0 directly on the full primitive/form closure."
                )
                if dq_zero_equivalent
                else (
                    "Characterize the primitive trace image Y=R({F: F'=original "
                    "test}) inside X_R and prove D_q|_Y=0, or prove direct "
                    "positivity of Q_Phi on the primitive image."
                )
            )
            if primitive_boundary_zero
            else (
                "Derive the exact primitive Weyl/Volterra boundary form "
                "beta_bdy=Q_original-Q_Phi, prove it descends to X_R, then "
                "compare the resulting D_bdy with D_q in transported trace "
                "coordinates."
            )
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Canonical boundary repair comparison")
    print(f"  D_q defined: {minimal_repair['closed']}")
    print(f"  canonical repair needed: {canonical_needed['closed']}")
    print(f"  abstract descent/comparison theorem: {abstract_descent['closed'] and comparison_sufficiency['closed']}")
    print(f"  Green identity input: {green_identity_input['closed']}")
    print(f"  D_bdy constructed: {exact_primitive_boundary['closed']}")
    print(f"  D_bdy descends through R: {boundary_descends['closed']}")
    print(f"  D_bdy >= D_q: {domination['closed']}")
    print(f"  original positivity closed: {original_positivity['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
