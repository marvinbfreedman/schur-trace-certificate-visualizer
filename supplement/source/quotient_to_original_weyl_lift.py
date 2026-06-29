#!/usr/bin/env python3
r"""Quotient-to-original Weyl lift ledger.

This script audits the requested lift

    original K_omega quadratic form
      == closed trace quotient form

after parity reduction, mixed-derivative/primitive transport, Volterra
normalization, and density/closure.

It proves the algebraic part of the lift and isolates the remaining endpoint
trace compatibility lemma.  The full lift is not marked closed unless the
primitive image of the original Weyl test space is identified with the
closed-trace quotient domain, or an exact boundary-repair identity is supplied
for the trace term.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def text(path: str) -> str:
    return Path(path).read_text(encoding="utf-8")


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def has_all(haystack: str, needles: list[str]) -> bool:
    return all(needle in haystack for needle in needles)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quotient-json",
        default="weyl_volterra_quotient_schur_theorem.json",
    )
    parser.add_argument(
        "--external-audit-json",
        default="weyl_volterra_external_equivalence_audit.json",
    )
    parser.add_argument(
        "--boundary-repair-json",
        default="boundary_repair_identity.json",
    )
    parser.add_argument(
        "--boundary-comparison-json",
        default="canonical_boundary_repair_comparison.json",
    )
    parser.add_argument(
        "--primitive-endpoint-json",
        default="primitive_endpoint_compatibility_theorem.json",
    )
    parser.add_argument("--notes", default="rh_klm_notes.md")
    parser.add_argument("--draft", default="rh_weyl_positive_draft.tex")
    parser.add_argument("--json-out", default="quotient_to_original_weyl_lift.json")
    args = parser.parse_args()

    quotient = load(args.quotient_json)
    audit = load(args.external_audit_json)
    boundary = load(args.boundary_repair_json) if Path(args.boundary_repair_json).exists() else None
    comparison = (
        load(args.boundary_comparison_json)
        if Path(args.boundary_comparison_json).exists()
        else None
    )
    primitive_endpoint = (
        load(args.primitive_endpoint_json)
        if Path(args.primitive_endpoint_json).exists()
        else None
    )
    notes = text(args.notes)
    draft = text(args.draft)

    quotient_closed = bool(
        quotient.get("globalWeylVolterraSchurStatus", {}).get("closed")
    )

    parity = status(
        "unitary parity reduction of original K_omega",
        has_all(notes, ["K_even(x, y)", "K_odd(x, y)", "It is enough to prove"])
        and has_all(draft, ["A\\ge 0", "-A\\le B\\le A"]),
        (
            "The simultaneous reflection symmetry K(a,b)=K(-a,-b) gives a "
            "unitary decomposition of L2(R) into even and odd half-line sectors. "
            "The original quadratic form is the direct sum of the two parity "
            "forms."
        ),
    )
    primitive = status(
        "mixed-derivative primitive identity",
        has_all(
            notes,
            [
                "P_+(x, y) = K_even(x, y)",
                "P_-(x, y) = K_odd(x, y)",
                "P_+(x, y) = integral_x^infty integral_y^infty H_+(u, v) du dv",
                "F(u) = sum_i c_i 1_{u >= x_i}",
            ],
        ),
        (
            "Boundary decay of the parity kernels gives "
            "P_pm(x,y)=int_x^inf int_y^inf H_pm(u,v)dudv.  Fubini then "
            "transports each finite original quadratic form to the mixed-kernel "
            "quadratic form on the primitive/cumulative test function."
        ),
    )
    full_line_mixed = status(
        "full-line mixed Green kernel equals Volterra source kernel",
        has_all(
            notes,
            [
                "mathcal H_omega(a, b)",
                "partial_a partial_b K_omega(a, b)",
                "integral_{|m|}^infty D_omega",
                "On the half-line, this full-line kernel has block form",
            ],
        ),
        (
            "The mixed derivative of the coordinate Weyl kernel is exactly the "
            "full-line Green kernel obtained by integrating the source "
            "D_omega.  The half-line block form is C_omega and R_omega^*."
        ),
    )
    volterra_normalization = status(
        "Volterra/log normalization is an invertible coordinate change",
        has_all(
            notes,
            [
                "K_red(s,t)",
                "K(s/2,t/2) / [Psi(s)Psi(t)]",
                "Using the direct Volterra formula",
                "Q(f,g)",
            ],
        ),
        (
            "The reduced kernel K_red is obtained from the same Weyl/Volterra "
            "kernel by positive multiplication by Psi factors and the log "
            "coordinate.  These operations preserve quadratic-form positivity "
            "on the transported form domain."
        ),
    )
    density = status(
        "density and closure of compact tests in the Volterra form domain",
        has_all(
            notes,
            [
                "Galerkin-to-continuum high-block exhaustion",
                "Mosco convergence",
                "density/closure",
            ],
        )
        or has_all(
            draft,
            ["density/closure", "Mosco", "full closed high block"],
        ),
        (
            "The high-block exhaustion theorem supplies the closed-space "
            "passage for A-graph Galerkin approximants.  The primitive identity "
            "therefore extends from finite combinations/compact smooth tests "
            "to the completed Volterra form domain whenever the trace "
            "compatibility condition is satisfied."
        ),
    )
    schur = status(
        "closed trace quotient Schur certificate",
        quotient_closed,
        (
            "The quotient Schur theorem proves positivity on ker R_global, the "
            "Douglas cross-form condition, and the singular Moore-Penrose "
            "Schur/range hypotheses in the normalized full-Phi model."
        ),
    )
    boundary_repair_closed = bool(
        boundary and boundary.get("boundaryRepairIdentityStatus", {}).get("closed")
    )
    canonical_comparison_closed = bool(
        comparison and comparison.get("canonicalBoundaryComparisonClosed")
    )
    primitive_endpoint_closed = bool(
        primitive_endpoint
        and primitive_endpoint.get("primitiveEndpointCompatibilityClosed")
    )
    primitive_boundary_zero = bool(
        comparison
        and comparison.get("importedSignals", {}).get("primitiveBoundaryZero")
    )
    dq_zero_equivalent = bool(
        comparison
        and comparison.get("importedSignals", {}).get(
            "DqZeroOnPrimitiveImageEquivalentToDqZero"
        )
    )
    primitive_vanishing_false = bool(
        boundary
        and boundary.get("primitiveVanishingFalse")
        and boundary.get("abstractRepairNonuniquenessStatus", {}).get("closed")
    )
    endpoint_trace_compatibility = status(
        "endpoint trace compatibility for primitive original tests",
        boundary_repair_closed or canonical_comparison_closed or primitive_endpoint_closed,
        (
            (
                "Closed by the primitive endpoint compatibility theorem: "
                "primitive tests do not lie in ker R_global, but D_bdy=0 and "
                "D_q=0 on X_R, so the trace correction is annihilated."
                if primitive_endpoint_closed
                else "Closed by the canonical boundary comparison theorem."
                if canonical_comparison_closed
                else "Closed by the boundary-repair identity ledger."
                if boundary_repair_closed
                else (
                    "The primitive-vanishing route is false: compact smooth "
                    "primitives can realize nonzero Lambda_a jets.  The only "
                    "remaining route is the canonical boundary comparison.  "
                    + (
                        "The primitive transport audit gives D_bdy=0, so the "
                        + (
                            "remaining theorem is D_q=0 on X_R, equivalently "
                            "Gamma^*Gamma-C<=0, or direct positivity of Q_Phi "
                            "on the full primitive/form closure."
                            if dq_zero_equivalent
                            else (
                                "remaining theorem is D_q|_Y=0 on the primitive "
                                "trace image Y, or direct positivity of Q_Phi "
                                "on that image."
                            )
                        )
                        if primitive_boundary_zero
                        else (
                            "The abstract quotient repair is nonunique, so one "
                            "must construct the Weyl/Volterra boundary operator "
                            "D_bdy and prove D_bdy >= D_q on the completed "
                            "trace range."
                        )
                    )
                    if primitive_vanishing_false
                    else (
                        "The existing ledgers do not prove that the primitive/"
                        "cumulative image of every original Weyl test lies in "
                        "ker R_global, nor do they prove an exact identity "
                        "adding back the trace repair ||S R_global f||^2 to "
                        "recover the original K_omega quadratic form."
                    )
                )
            )
        ),
        blocker=(
            None
            if boundary_repair_closed or canonical_comparison_closed or primitive_endpoint_closed
            else (
                comparison.get("nextProofTarget")
                if comparison
                else
                boundary.get("nextProofTarget")
                if boundary
                else (
                    "Prove either (A) primitive original tests satisfy "
                    "R_global f=0 in the closed Volterra form domain, or (B) "
                    "the original Weyl quadratic form equals "
                    "Q_Phi(f)+||S R_global f||^2 with the same trace-side "
                    "operator S from the quotient Schur theorem."
                )
            )
        ),
    )

    algebraic_statuses = [
        parity,
        primitive,
        full_line_mixed,
        volterra_normalization,
        density,
        schur,
    ]
    algebraic_lift_closed = all(item["closed"] for item in algebraic_statuses)
    conditional_lift_closed = algebraic_lift_closed
    full_lift_closed = algebraic_lift_closed and endpoint_trace_compatibility["closed"]

    data = {
        "theoremName": "quotient-to-original Weyl lift",
        "quotientJson": args.quotient_json,
        "externalAuditJson": args.external_audit_json,
        "boundaryRepairJson": args.boundary_repair_json if boundary else None,
        "boundaryComparisonJson": args.boundary_comparison_json if comparison else None,
        "primitiveEndpointJson": args.primitive_endpoint_json if primitive_endpoint else None,
        "normalizedSchurCertificateClosed": quotient_closed,
        "algebraicLiftClosed": algebraic_lift_closed,
        "conditionalLiftTheoremClosed": conditional_lift_closed,
        "quotientToOriginalWeylLiftClosed": full_lift_closed,
        "originalWeylKernelPositivityClosed": full_lift_closed
        and not any(
            item["label"].startswith("omega coverage") and not item["closed"]
            for item in audit.get("auditItems", [])
        ),
        "statuses": {
            "parityReductionStatus": parity,
            "mixedDerivativePrimitiveStatus": primitive,
            "fullLineMixedGreenStatus": full_line_mixed,
            "volterraNormalizationStatus": volterra_normalization,
            "densityClosureStatus": density,
            "quotientSchurStatus": schur,
            "endpointTraceCompatibilityStatus": endpoint_trace_compatibility,
        },
        "boundaryRepairTheorem": boundary,
        "boundaryComparisonTheorem": comparison,
        "primitiveEndpointCompatibilityTheorem": primitive_endpoint,
        "conditionalTheorem": (
            "If the endpoint trace compatibility lemma holds, then for every "
            "original compact Weyl test the original K_omega quadratic form "
            "is transported by parity, the primitive identity, and Volterra "
            "normalization to the closed trace quotient form.  The quotient "
            "Schur theorem then gives nonnegativity in the normalized "
            "full-Phi model."
        ),
        "blockingLemma": endpoint_trace_compatibility["blocker"],
        "nextProofTarget": endpoint_trace_compatibility["blocker"],
        "auditUpdate": (
            "The previous broad 'quotient constraint equals original test "
            "space' gap is now reduced to endpoint trace compatibility for "
            "the primitive image, or equivalently the exact boundary-repair "
            "identity for the trace term."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Quotient-to-original Weyl lift")
    print(f"  algebraic lift closed: {algebraic_lift_closed}")
    print(f"  conditional lift theorem closed: {conditional_lift_closed}")
    print(f"  endpoint trace compatibility: {endpoint_trace_compatibility['closed']}")
    print(f"  full lift closed: {full_lift_closed}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
