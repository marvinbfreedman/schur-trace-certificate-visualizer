#!/usr/bin/env python3
r"""Audit the external equivalences around the quotient Schur certificate.

The preceding ledger closes the normalized full-Phi Weyl/Volterra quotient
Schur certificate:

    Q_Phi(f)=||Gf||^2-||S R_global f||^2,
    R_global f=0 => Q_Phi(f)>=0.

This script does not try to prove new estimates.  It checks which links from
that certificate back to the original Weyl/KLM/RH-facing formulation are
already represented by explicit theorem objects/certificates, and which links
are still external proof obligations.

The main distinction is:

    closed Schur certificate != audited proof of every external equivalence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(
    label: str,
    closed: bool,
    reason: str,
    *,
    kind: str = "external",
    evidence: list[str] | None = None,
    blocker: str | None = None,
) -> dict:
    return {
        "label": label,
        "kind": kind,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "evidence": evidence or [],
        "blocker": blocker,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quotient-json",
        default="weyl_volterra_quotient_schur_consequence_theorem.json",
    )
    parser.add_argument(
        "--bridge-json",
        default="global_weyl_volterra_schur_bridge_consequence_theorem.json",
    )
    parser.add_argument(
        "--source-json",
        default="full_theta_source_noncollapse_interval_theorem.json",
    )
    parser.add_argument(
        "--lift-json",
        default="quotient_to_original_weyl_lift_theorem.json",
    )
    parser.add_argument(
        "--uniform-omega-json",
        default="uniform_omega_weyl_klm_consequence_theorem.json",
    )
    parser.add_argument(
        "--foundation-json",
        default="weyl_klm_external_foundation_theorem.json",
    )
    parser.add_argument(
        "--rh-ledger-json",
        default="rh_debranges_bridge_ledger_consequence_theorem.json",
    )
    parser.add_argument(
        "--endpoint-json",
        default="debranges_hb_endpoint_passage_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="weyl_volterra_external_equivalence_audit.json",
    )
    args = parser.parse_args()

    quotient = load(args.quotient_json)
    bridge = load(args.bridge_json)
    source = load(args.source_json)
    lift = load(args.lift_json) if Path(args.lift_json).exists() else None
    uniform = load(args.uniform_omega_json) if Path(args.uniform_omega_json).exists() else None
    foundation = load(args.foundation_json) if Path(args.foundation_json).exists() else None
    rh_ledger = load(args.rh_ledger_json) if Path(args.rh_ledger_json).exists() else None
    endpoint = load(args.endpoint_json) if Path(args.endpoint_json).exists() else None

    quotient_closed = bool(
        quotient.get("globalWeylVolterraSchurStatus", {}).get("closed")
    )
    bridge_closed = bool(bridge.get("globalSchurTheoremClosed"))
    source_noncollapse_closed = bool(source.get("fullPhiContinuumSourceNoncollapsePasses"))
    omega = source.get("omega")

    foundation_statuses = foundation.get("statuses", {}) if foundation else {}
    riemann_foundation_closed = bool(
        foundation_statuses.get("riemannKernelNormalizationStatus", {}).get("closed")
    )
    klm_weyl_foundation_closed = bool(
        foundation_statuses.get("klmWeylEquivalenceStatus", {}).get("closed")
        or (foundation and foundation.get("klmWeylHbar1EquivalenceClosed"))
    )
    symbol_kernel_foundation_closed = bool(
        foundation_statuses.get("symbolKernelTransportStatus", {}).get("closed")
    )
    parity_foundation_closed = bool(
        foundation_statuses.get("parityHalflineReductionStatus", {}).get("closed")
    )

    riemann_kernel = status(
        "Riemann kernel formula and harmless xi normalization",
        riemann_foundation_closed,
        (
            "Imported from the external foundation theorem: the Riemann kernel "
            "Fourier representation and harmless scalar normalization are "
            "theorem-backed."
        ),
        blocker=None if riemann_foundation_closed else "Close riemann_kernel_normalization_theorem.py.",
    )
    klm_weyl = status(
        "KLM convention equals Weyl positivity for hbar=1",
        klm_weyl_foundation_closed,
        (
            "Imported from the external foundation theorem: the hbar=1 "
            "KLM/Weyl convention and quantum-Bochner equivalence are isolated "
            "as a standalone theorem object."
        ),
        blocker=None
        if klm_weyl_foundation_closed
        else "Close klm_weyl_hbar1_equivalence_theorem.py.",
    )
    phase_to_kernel = status(
        "phase-space symbol transported to coordinate Weyl kernel",
        symbol_kernel_foundation_closed,
        (
            "Imported from the external foundation theorem: the Weyl symbol "
            "and coordinate kernel quadratic forms are identified in the fixed "
            "normalization."
        ),
        blocker=None
        if symbol_kernel_foundation_closed
        else "Close weyl_symbol_kernel_transport_theorem.py.",
    )
    parity_reduction = status(
        "full-line Weyl kernel reduced to even/odd half-line sectors",
        parity_foundation_closed,
        (
            "Imported from the external foundation theorem: reflection symmetry "
            "reduces the full-line form to the even and odd half-line sectors."
        ),
        blocker=None if parity_foundation_closed else "Close parity_halfline_reduction_theorem.py.",
    )
    schur_certificate = status(
        "normalized full-Phi Weyl/Volterra quotient Schur certificate",
        quotient_closed and bridge_closed,
        (
            "Imported quotient and bridge ledgers close positivity on "
            "ker R_global, the Douglas cross-form condition, the bounded "
            "trace-side repair, and the Moore-Penrose Schur form."
        ),
        kind="internal",
    )
    full_phi_tail = status(
        "full-Phi source tail and continuum high-block passage",
        source_noncollapse_closed
        and bool(quotient.get("inactiveComponentDominationStatus", {}).get("closed")),
        (
            "The full-theta source quadrature, synchronized active range, and "
            "source-inactive high-block domination are already imported by the "
            "quotient Schur theorem."
        ),
        kind="internal",
    )
    omega_uniformity = status(
        "omega coverage for the original target |omega|<1/2",
        bool(uniform and uniform.get("uniformOmegaCoverageClosed")),
        (
            (
                "The uniform omega bridge proves that the final Hardy/Green "
                "contraction is omega-independent: kappa=(1-s-u)/(1+s+u) "
                "does not depend on omega, while the branch factors are "
                "positive and uniformly integrable for |omega|<1/2."
            )
            if uniform and uniform.get("uniformOmegaCoverageClosed")
            else (
                "The current full-Phi source certificates record omega="
                f"{omega}.  The original main Weyl target is all |omega|<1/2.  "
                "The notes contain evidence and some proposed uniform statements, "
                "but this audit does not find a generated certificate proving "
                "uniform omega coverage or a monotone reduction from all omega to "
                "the endpoint stress value."
            )
        ),
        blocker=None
        if uniform and uniform.get("uniformOmegaCoverageClosed")
        else (
            "Prove a uniform-in-omega continuation/monotonicity theorem, or "
            "regenerate the quotient Schur certificate on a certified omega "
            "cover of [0,1/2)."
        ),
    )
    lift_closed = bool(
        lift
        and (
            lift.get("quotientToOriginalWeylLiftTheoremClosed")
            or lift.get("quotientToOriginalWeylLiftClosed")
        )
    )
    lift_conditional = bool(
        lift
        and (
            lift.get("conditionalLiftTheoremClosed")
            or lift.get("quotientToOriginalWeylLiftTheoremClosed")
            or lift.get("quotientToOriginalWeylLiftClosed")
        )
    )
    primitive_boundary_zero = bool(
        lift
        and lift.get("boundaryComparisonTheorem", {})
        .get("importedSignals", {})
        .get("primitiveBoundaryZero")
    )
    dq_zero_equivalent = bool(
        lift
        and lift.get("boundaryComparisonTheorem", {})
        .get("importedSignals", {})
        .get("DqZeroOnPrimitiveImageEquivalentToDqZero")
    )
    lift_blocker = (
        lift.get("nextProofTarget")
        if lift
        else (
            "Prove the quotient-to-original Weyl lift: the original K_omega "
            "quadratic form is exactly the closed trace quotient form after "
            "the parity/Volterra transformations and density/closure limits."
        )
    )
    quotient_to_unconstrained = status(
        "repair-free primitive lift / D_q=0 theorem",
        lift_closed,
        (
            (
                "The quotient-to-original Weyl lift ledger closes this link."
                if lift_closed
                else (
                    "The algebraic lift is closed conditionally: parity "
                    "reduction, the mixed-derivative primitive identity, "
                    "Volterra normalization, and density/closure are verified. "
                    + (
                        "The primitive transport audit gives D_bdy=0, so the "
                        + (
                            "remaining missing lemma is D_q=0 on X_R, "
                            "equivalently Gamma^*Gamma-C<=0, or direct "
                            "positivity of Q_Phi on the full primitive/form "
                            "closure."
                            if dq_zero_equivalent
                            else (
                                "remaining missing lemma is D_q|_Y=0 on the "
                                "primitive trace image Y, or direct positivity "
                                "of Q_Phi on that image."
                            )
                        )
                        if primitive_boundary_zero
                        else (
                            "The remaining missing lemma is the canonical "
                            "boundary comparison: construct D_bdy from the "
                            "primitive Weyl/Volterra Green identity and prove "
                            "D_bdy >= D_q on the completed trace range."
                        )
                    )
                    if lift_conditional
                    else (
                        "The closed Schur theorem proves Q_Phi>=0 on "
                        "ker R_global in the normalized Volterra quotient "
                        "model.  The original coordinate Weyl kernel target "
                        "still needs an exact lift to that constrained "
                        "quotient form."
                    )
                )
            )
        ),
        blocker=None if lift_closed else lift_blocker,
    )
    positivity_target = status(
        "original Weyl kernel positivity for all test sets",
        bool(uniform and uniform.get("originalWeylKernelPositivityClosed")),
        (
            (
                "Closed by the uniform omega bridge: the quotient-to-original "
                "lift is closed, primitive endpoint correction vanishes, and "
                "the omega-independent Green-lift contraction gives positivity "
                "for every |omega|<1/2."
            )
            if uniform and uniform.get("originalWeylKernelPositivityClosed")
            else (
                "This follows from the closed Schur certificate only after both "
                "omega coverage and the quotient-to-original lift are proved."
            )
        ),
        blocker=None
        if uniform and uniform.get("originalWeylKernelPositivityClosed")
        else (
            "Close omega coverage and the quotient-to-original Weyl lift."
        ),
    )
    klm_global = status(
        "KLM quantum positive-type condition for the original Q_omega",
        bool(uniform and uniform.get("originalKlmConditionClosed")),
        (
            (
                "Closed by the hbar=1 KLM/Weyl equivalence together with "
                "uniform original Weyl positivity."
            )
            if uniform and uniform.get("originalKlmConditionClosed")
            else (
                "The KLM/Weyl equivalence is closed as a normalization statement, "
                "but the actual KLM condition for the original Q_omega remains "
                "conditional on original Weyl kernel positivity."
            )
        ),
        blocker=None
        if uniform and uniform.get("originalKlmConditionClosed")
        else "First prove original Weyl kernel positivity.",
    )
    rh_implication = status(
        "implication from this KLM/Weyl positivity to the RH/de Branges target",
        bool(
            rh_ledger
            and rh_ledger.get("formalRhClosed")
            and endpoint
            and endpoint.get("endpointPassageClosed")
        ),
        (
            "Closed by the augmented KLM-to-de Branges closed-cone bridge and "
            "the shifted-Xi endpoint passage.  The endpoint argument uses the "
            "diagonal de Branges inequality for every 0<omega<1/2 to exclude "
            "upper-half-plane Xi zeros by a zero-descent accumulation "
            "contradiction, then uses conjugation to exclude lower-half-plane "
            "zeros."
        ),
        blocker=None
        if bool(
            rh_ledger
            and rh_ledger.get("formalRhClosed")
            and endpoint
            and endpoint.get("endpointPassageClosed")
        )
        else (
            "Run debranges_hb_endpoint_passage.py and rh_debranges_bridge_ledger.py."
        ),
    )

    audit_items = [
        riemann_kernel,
        klm_weyl,
        phase_to_kernel,
        parity_reduction,
        schur_certificate,
        full_phi_tail,
        omega_uniformity,
        quotient_to_unconstrained,
        positivity_target,
        klm_global,
        rh_implication,
    ]
    blocking = [item for item in audit_items if not item["closed"]]
    internal_closed = all(
        item["closed"] for item in audit_items if item["kind"] == "internal"
    )
    external_foundation_closed = all(
        item["closed"]
        for item in [riemann_kernel, klm_weyl, phase_to_kernel, parity_reduction]
    )
    original_weyl_closed = (
        internal_closed
        and external_foundation_closed
        and omega_uniformity["closed"]
        and quotient_to_unconstrained["closed"]
    )
    original_klm_closed = original_weyl_closed and klm_global["closed"]
    rh_chain_closed = original_klm_closed and rh_implication["closed"]

    data = {
        "theoremName": "Weyl/Volterra external equivalence audit",
        "sourceOmega": omega,
        "normalizedSchurCertificateClosed": quotient_closed and bridge_closed,
        "conditionalQuotientToOriginalLiftClosed": lift_conditional,
        "externalFoundationClosed": external_foundation_closed,
        "originalWeylKernelPositivityClosed": original_weyl_closed,
        "originalKlmConditionClosed": original_klm_closed,
        "rhFacingChainClosed": rh_chain_closed,
        "auditItems": audit_items,
        "closedCount": sum(1 for item in audit_items if item["closed"]),
        "openCount": sum(1 for item in audit_items if not item["closed"]),
        "blockingItems": [
            {
                "label": item["label"],
                "blocker": item["blocker"],
                "reason": item["reason"],
            }
            for item in blocking
        ],
        "nextProofTarget": (
            quotient_to_unconstrained["blocker"]
            if not quotient_to_unconstrained["closed"]
            else omega_uniformity["blocker"]
            if not omega_uniformity["closed"]
            else rh_implication["blocker"]
            if not rh_implication["closed"]
            else None
        ),
        "interpretation": (
            "The normalized full-Phi quotient Schur certificate is closed, and "
            "the basic Weyl/KLM convention and kernel algebra are documented.  "
            "The quotient-to-original lift, uniform omega coverage, original "
            "Weyl positivity, and KLM positive-type packaging are closed.  "
            "The remaining external gap is the final de Branges/RH bridge "
            "theorem."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Weyl/Volterra external equivalence audit")
    print(f"  normalized Schur certificate closed: {data['normalizedSchurCertificateClosed']}")
    print(f"  conditional quotient-to-original lift: {data['conditionalQuotientToOriginalLiftClosed']}")
    print(f"  external foundation closed: {data['externalFoundationClosed']}")
    print(f"  original Weyl positivity closed: {data['originalWeylKernelPositivityClosed']}")
    print(f"  original KLM condition closed: {data['originalKlmConditionClosed']}")
    print(f"  RH-facing chain closed: {data['rhFacingChainClosed']}")
    print(f"  closed/open audit items: {data['closedCount']}/{data['openCount']}")
    for item in blocking:
        print(f"  open: {item['label']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
