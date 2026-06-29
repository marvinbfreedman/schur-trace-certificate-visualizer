#!/usr/bin/env python3
r"""Primitive endpoint compatibility via the Green-lift contraction.

The original Weyl primitive route had two possible endpoints:

  (A) primitive original tests land in ker R_global;
  (B) the boundary/trace correction from transport dominates the quotient
      repair D_q.

The existing primitive ledgers prove that (A) is false and that the primitive
transport boundary operator is actually zero:

    D_bdy = 0.

Therefore (B) reduces to the sharper repair-free statement

    D_q = 0 on X_R.

This theorem closes that statement on the completed Volterra trace-fiber
domain by identifying the quotient Schur trace form with the signed
Volterra/Green feature form

    D_trace = P - M,

where P is the plus-feature Gram and M is the minus-feature Gram on Green
minimizer trace images.  The continuum Green-lift closure theorem proves

    M <= P

because G_- = C K E G_+ and ||C K E|| <= 1.  Hence the trace-side Schur form
P-M is nonnegative, the positive part (M-P)_+ vanishes, and therefore

    D_q = 0.

Since the primitive trace image is dense in X_R and D_q is bounded, this also
proves D_q vanishes on primitive traces.  Together with D_bdy=0, the primitive
endpoint compatibility layer is closed in the completed Volterra model.
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
    parser.add_argument("--primitive-boundary-json", default="primitive_boundary_transport_audit.json")
    parser.add_argument("--primitive-density-json", default="primitive_trace_image_density.json")
    parser.add_argument("--feature-json", default="trace_volterra_green_feature_map.json")
    parser.add_argument("--closure-json", default="continuum_green_lift_closure_theorem.json")
    parser.add_argument(
        "--quotient-json",
        default="quotient_minimal_repair_consequence_theorem.json",
    )
    parser.add_argument(
        "--consequence-json",
        default="primitive_endpoint_compatibility_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="primitive_endpoint_compatibility_theorem.json")
    args = parser.parse_args()

    primitive_boundary = load(args.primitive_boundary_json)
    primitive_density = load(args.primitive_density_json)
    feature = load(args.feature_json)
    closure = load(args.closure_json)
    quotient = load(args.quotient_json)
    consequence = load(args.consequence_json)
    consequence_closed = bool(consequence.get("primitiveEndpointCompatibilityConsequenceClosed"))

    primitive_a_false = nested_closed(
        primitive_boundary, "statuses", "comparisonReducesToDqZeroStatus"
    )
    dbdy_zero = nested_closed(
        primitive_boundary, "statuses", "canonicalPrimitiveBoundaryZeroStatus"
    )
    dbdy_descends = nested_closed(
        primitive_boundary, "statuses", "zeroBoundaryDescendsStatus"
    )
    primitive_dense = nested_closed(
        primitive_density, "statuses", "primitiveTraceDenseStatus"
    )
    dq_equiv = nested_closed(
        primitive_density, "statuses", "dqZeroEquivalenceStatus"
    )
    volterra_identity = nested_closed(
        feature, "statuses", "volterraMomentRepresentationStatus"
    )
    signed_square = nested_closed(
        feature, "statuses", "signedSquareCompletionStatus"
    )
    green_contraction = nested_closed(
        closure, "statuses", "greenLiftContractionStatus"
    )
    quotient_repair_bounded = bool(
        quotient.get("traceSideRepairClosed")
        or nested_closed(quotient, "traceSideRepairStatus")
    )
    quotient_factorization = bool(
        quotient.get("quotientFactorizationClosed")
        or nested_closed(quotient, "quotientFactorizationStatus")
    )

    route_a_status = status(
        "primitive tests land in ker R_global",
        False,
        (
            "False.  The primitive boundary ledger gives a compact smooth "
            "jet-extension counterexample with Lambda_a(F) != 0 for an active "
            "trace row."
        ),
        blocker="Use the repair-free route instead.",
    )
    zero_boundary_status = status(
        "primitive transport boundary repair is zero",
        dbdy_zero and dbdy_descends,
        (
            "The half-line primitive integration-by-parts identity has no "
            "endpoint contribution: F(0)=0 kills the lower endpoint and kernel "
            "decay kills infinity.  Thus beta_bdy=0 and D_bdy=0; the zero "
            "form descends through R_global."
        ),
        blocker=None if dbdy_zero and dbdy_descends else "Close primitive_boundary_transport_audit.py.",
    )
    feature_identification_status = status(
        "quotient trace Schur form equals plus/minus Green feature form",
        volterra_identity and signed_square and quotient_factorization,
        (
            "On Green-minimizer trace images, the transported trace form is "
            "D_trace=P-M, where P=<G_+,G_+> and M=<G_-,G_->.  This is the "
            "same Schur trace form used by the quotient factorization."
        ),
        blocker=None
        if volterra_identity and signed_square and quotient_factorization
        else "Close the Volterra Green feature-map identity and quotient factorization.",
    )
    minus_dominated_status = status(
        "minus feature Gram dominated by plus feature Gram",
        green_contraction,
        (
            "The continuum Green-lift closure theorem gives "
            "G_-=C K E G_+ with ||C K E||<=1 on the completed Green-minimizer "
            "trace image.  Hence M=<G_-,G_-> <= <G_+,G_+>=P."
        ),
        blocker=None if green_contraction else "Close continuum_green_lift_closure_theorem.py.",
    )
    dq_zero_status = status(
        "quotient repair D_q vanishes on X_R",
        feature_identification_status["closed"]
        and minus_dominated_status["closed"]
        and quotient_repair_bounded,
        (
            "Since D_trace=P-M>=0, the positive repair part "
            "D_q=(M-P)_+ vanishes in the plus/minus Green feature coordinates. "
            "Equivalently Gamma^*Gamma<=C on X_R."
        ),
        blocker=None
        if feature_identification_status["closed"] and minus_dominated_status["closed"] and quotient_repair_bounded
        else "Need feature identification, Green-lift contraction, and bounded quotient repair.",
    )
    primitive_trace_status = status(
        "D_q vanishes on primitive trace image",
        dq_zero_status["closed"] and primitive_dense and dq_equiv,
        (
            "The primitive trace image is dense in X_R and D_q is bounded.  "
            "Thus D_q=0 on X_R is equivalent to D_q=0 on primitive traces."
        ),
        blocker=None
        if dq_zero_status["closed"] and primitive_dense and dq_equiv
        else "Close primitive trace density and D_q=0 on X_R.",
    )
    endpoint_compat_status = status(
        "primitive endpoint compatibility",
        zero_boundary_status["closed"] and primitive_trace_status["closed"],
        (
            "Although primitive tests do not lie in ker R_global, the boundary "
            "repair from primitive transport is zero and the quotient repair "
            "also vanishes on the completed trace range.  Therefore the trace "
            "correction is annihilated on primitive traces."
        ),
        blocker=None
        if zero_boundary_status["closed"] and primitive_trace_status["closed"]
        else "Close D_bdy=0 and D_q=0.",
    )

    data = {
        "theoremName": "primitive endpoint compatibility theorem",
        "proofClass": "symbolic identity",
        "consequenceJson": args.consequence_json,
        "primitiveBoundaryJson": args.primitive_boundary_json,
        "primitiveDensityJson": args.primitive_density_json,
        "featureJson": args.feature_json,
        "closureJson": args.closure_json,
        "quotientJson": args.quotient_json,
        "importedSignals": {
            "primitiveRouteAFalse": primitive_a_false,
            "D_bdyZero": dbdy_zero,
            "D_bdyDescends": dbdy_descends,
            "primitiveTraceDense": primitive_dense,
            "DqZeroOnPrimitiveImageEquivalentToDqZero": dq_equiv,
            "volterraMomentRepresentation": volterra_identity,
            "signedSquareCompletion": signed_square,
            "greenLiftContraction": green_contraction,
            "quotientRepairBounded": quotient_repair_bounded,
            "quotientFactorization": quotient_factorization,
        },
        "statuses": {
            "routeAStatus": route_a_status,
            "zeroPrimitiveBoundaryRepairStatus": zero_boundary_status,
            "greenFeatureSchurIdentificationStatus": feature_identification_status,
            "minusDominatedByPlusStatus": minus_dominated_status,
            "dqVanishesOnXRStatus": dq_zero_status,
            "dqVanishesOnPrimitiveTraceStatus": primitive_trace_status,
            "primitiveEndpointCompatibilityStatus": endpoint_compat_status,
        },
        "operatorIdentifications": {
            "D_bdy": "0",
            "D_trace": "P-M on Green-minimizer trace images",
            "P": "<G_+,G_+>",
            "M": "<G_-,G_->",
            "contraction": "G_- = C K E G_+, ||C K E|| <= 1",
            "D_q": "(M-P)_+ = 0, equivalently Gamma^*Gamma-C <= 0",
        },
        "proof": [
            "Primitive route A is false: compact primitives can have nonzero active trace.",
            "Primitive transport has no endpoint boundary term, so D_bdy=0.",
            "The primitive trace image is dense in X_R, so D_q|Y=0 iff D_q=0 on X_R.",
            "On Green-minimizer trace images, the Schur trace form is P-M.",
            "The Green-lift closure theorem gives M<=P.",
            "Therefore D_q=(M-P)_+=0 on X_R and hence on primitive traces.",
        ],
        "primitiveEndpointCompatibilityClosed": endpoint_compat_status["closed"],
        "primitiveEndpointCompatibilityConsequenceClosed": consequence_closed,
        "formalConclusion": (
            "Original Weyl primitives do not generally land in ker R_global.  "
            "However, the primitive boundary correction is zero and the "
            "quotient repair vanishes because the completed Green-lift "
            "contraction gives M<=P.  Thus the primitive endpoint trace "
            "correction is annihilated in the completed Volterra model."
        ),
        "nextProofTarget": (
            "Update the quotient-to-original Weyl lift ledger to import this "
            "primitive endpoint compatibility theorem."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive endpoint compatibility theorem")
    print(f"  route A lands in ker R: {route_a_status['closed']} (false route)")
    print(f"  D_bdy=0: {zero_boundary_status['closed']}")
    print(f"  feature Schur identification: {feature_identification_status['closed']}")
    print(f"  M <= P by Green-lift contraction: {minus_dominated_status['closed']}")
    print(f"  D_q=0 on X_R: {dq_zero_status['closed']}")
    print(f"  primitive endpoint compatibility: {endpoint_compat_status['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
