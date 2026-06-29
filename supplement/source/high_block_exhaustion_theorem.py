#!/usr/bin/env python3
r"""High-block Galerkin-to-continuum exhaustion theorem ledger.

This records the proof layer needed after the normalized source-inactive
tail constants consequence.  The certified source model proves

    ||(I-P_2) Ehat_N f||^2 <= eps_delta <A_N f,f>.

To upgrade this to the continuum closed-trace high block, one needs the
Mosco/strong-resolvent passage

    H_{M,N} -> H_M cap ker R_global

in the A graph norm and norm convergence of the compact source operators

    S_N = Ehat_N^* Ehat_N -> S = Ehat_Phi^* Ehat_Phi.

The script deliberately separates the formal functional-analytic implication
from the remaining analytic elliptic/Mosco hypotheses.  It must not mark the
global Schur theorem closed unless those hypotheses are discharged.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.is_file():
        return None
    return json.loads(candidate.read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def theorem_summary(data: dict | None, keys: list[str]) -> dict:
    if not data:
        return {}
    out = {"theoremName": data.get("theoremName")}
    for key in keys:
        if key in data:
            out[key] = data[key]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tail-constants-json",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument("--minmax-json", default="")
    parser.add_argument(
        "--elliptic-json",
        default="",
    )
    parser.add_argument(
        "--compact-consequence-json",
        default="high_block_compact_exhaustion_consequence_theorem.json",
    )
    parser.add_argument("--compact-proof-json", default="")
    parser.add_argument("--frame-json", default="")
    parser.add_argument("--quadrature-json", default="")
    parser.add_argument("--json-out", default="high_block_exhaustion_theorem.json")
    args = parser.parse_args()

    tail_constants_path = args.tail_constants_json or args.minmax_json
    minmax = load_optional(tail_constants_path)
    if not minmax:
        raise SystemExit(f"missing required tail constants certificate: {tail_constants_path}")

    elliptic_data = load_optional(args.elliptic_json)
    compact_consequence = load_optional(args.compact_consequence_json)
    compact_proof = load_optional(args.compact_proof_json)
    frame = load_optional(args.frame_json)
    quadrature = load_optional(args.quadrature_json)

    compact_conditional_closed = bool(
        (
            compact_consequence
            and compact_consequence.get("conditionalHighBlockExhaustionClosed")
        )
        or (compact_proof and compact_proof.get("conditionalHighBlockExhaustionClosed"))
    )
    compact_continuum_closed = bool(
        (
            compact_consequence
            and compact_consequence.get("tailEstimatePassesToContinuum")
        )
        or (compact_proof and compact_proof.get("tailEstimatePassesToContinuum"))
    )
    elliptic_closed = bool(
        compact_continuum_closed
        or (elliptic_data and elliptic_data.get("ellipticTraceEstimateClosed"))
    )
    continuum_frame_status = (
        (compact_proof or {}).get("continuumTraceFrameLowerBoundStatus", {})
        if compact_proof
        else {}
    )
    compact_reason = (
        (compact_proof or {}).get("remainingAnalyticGap")
        if compact_proof
        else None
    )
    elliptic_reason = (
        "The compact-source exhaustion proof supplies the valid replacement "
        "for the invalid raw compact-kernel Sobolev route."
        if compact_continuum_closed
        else (
            compact_reason
            or elliptic_data.get("remainingAnalyticGap")
            if elliptic_data
            else (
                "Needed to control local H^m jets from the A graph norm on "
                "H_M cap ker R_global; no elliptic trace ledger was found."
            )
        )
    )
    elliptic = status(
        "commuted Sturm/Volterra elliptic estimate",
        elliptic_closed,
        elliptic_reason,
    )
    mosco_limsup = status(
        "Mosco limsup/recovery sequence",
        compact_continuum_closed,
        (
            "The compact-source proof gives the recovery sequence once the "
            "sampled trace map has a continuum uniform lower frame bound and "
            "bounded correction right inverse."
        ),
    )
    mosco_liminf = status(
        "Mosco liminf/compactness",
        compact_continuum_closed,
        (
            "The compact-source proof gives the weak compactness/trace-closure "
            "argument once sampled trace quadrature is known to converge on "
            "A-bounded high-block sequences."
        ),
    )
    source_norm = status(
        "source operator norm convergence",
        compact_continuum_closed,
        (
            "The compact-source proof shows that Mosco convergence of the "
            "closed high blocks plus compactness of S=Ehat^*Ehat gives "
            "||Pi_N S Pi_N-S||_{A->A}->0."
        ),
    )
    abstract_compact_exhaustion = status(
        "abstract compact-source exhaustion theorem",
        compact_conditional_closed,
        (
            "The Hilbert-space Mosco/compact-source/min-max implication is "
            "closed.  It is conditional only on the continuum trace-frame lower "
            "bound needed to identify the sampled closed-trace Galerkin spaces "
            "with H_M."
        ),
    )
    spectral_projection = status(
        "Riesz spectral projection convergence",
        True,
        (
            "Standard consequence of ||P_N S P_N-S||_{A->A}->0 and an open "
            "gap around the top-two source spectrum."
        ),
    )
    minmax_passage = status(
        "min-max inequality passage",
        True,
        (
            "Once source operator norm convergence and spectral projection "
            "convergence hold, the Courant-Fischer/Riesz projection argument "
            "passes the normalized tail inequality to the continuum."
        ),
    )

    analytic_statuses = [
        elliptic,
        mosco_limsup,
        mosco_liminf,
        source_norm,
        spectral_projection,
        minmax_passage,
    ]
    analytic_hypotheses_closed = all(item["closed"] for item in analytic_statuses[:4])
    conditional_passage_closed = spectral_projection["closed"] and minmax_passage["closed"]
    tail_passes = analytic_hypotheses_closed and conditional_passage_closed

    finite_evidence = {}
    if frame:
        finite_evidence["observedGammaFloor"] = frame.get("observedGammaFloor")
        finite_evidence["maxRangeResidualRelative"] = frame.get("maxRangeResidualRelative")
        finite_evidence["maxDensityOperator"] = frame.get("maxDensityOperator")
        finite_evidence["maxRowDensityL2"] = frame.get("maxRowDensityL2")
        finite_evidence["frameInputs"] = frame.get("inputs")
    if quadrature:
        finite_evidence["maxAbsFrameMinRelativeDrift"] = quadrature.get(
            "maxAbsFrameMinRelativeDrift"
        )
        finite_evidence["maxAbsMaxRowDensityL2RelativeDrift"] = quadrature.get(
            "maxAbsMaxRowDensityL2RelativeDrift"
        )

    data = {
        "theoremName": "Galerkin-to-continuum high-block exhaustion",
        "tailConstantsJson": tail_constants_path,
        "legacyMinmaxJson": args.minmax_json or None,
        "diagnosticEllipticJson": args.elliptic_json if elliptic_data else None,
        "compactConsequenceJson": args.compact_consequence_json
        if compact_consequence
        else None,
        "legacyCompactProofJson": args.compact_proof_json if compact_proof else None,
        "frameJson": args.frame_json if frame else None,
        "quadratureJson": args.quadrature_json if quadrature else None,
        "normalizedEpsilonDelta": minmax["normalizedEpsilonDelta"],
        "finiteLowMidSchurBudget": minmax["finiteLowMidSchurBudget"],
        "absorptionSlack": minmax["absorptionSlack"],
        "sourceTopLower": minmax["sourceTopLower"],
        "absoluteTailBound": minmax["absoluteTailBound"],
        "ellipticTraceTheoremSummary": theorem_summary(
            elliptic_data,
            [
                "ellipticTraceEstimateClosed",
                "highBlockExhaustionInputClosed",
                "remainingAnalyticGap",
            ],
        ),
        "compactSourceExhaustionConsequenceSummary": theorem_summary(
            compact_consequence,
            [
                "conditionalHighBlockExhaustionClosed",
                "tailEstimatePassesToContinuum",
                "remainingAnalyticGap",
            ],
        ),
        "legacyCompactSourceExhaustionProofSummary": theorem_summary(
            compact_proof,
            [
                "conditionalHighBlockExhaustionClosed",
                "tailEstimatePassesToContinuum",
                "remainingAnalyticGap",
            ],
        ),
        "abstractCompactExhaustionStatus": abstract_compact_exhaustion,
        "continuumTraceFrameLowerBoundStatus": continuum_frame_status,
        "moscoLimsupStatus": mosco_limsup,
        "moscoLiminfStatus": mosco_liminf,
        "ellipticTraceStatus": elliptic,
        "sourceOperatorNormConvergenceStatus": source_norm,
        "spectralProjectionConvergenceStatus": spectral_projection,
        "minmaxPassageStatus": minmax_passage,
        "conditionalTailEstimatePassage": conditional_passage_closed,
        "tailEstimatePassesToContinuum": tail_passes,
        "globalSchurTheoremClosed": False,
        "finiteEvidence": finite_evidence,
        "formalProof": [
            (
                "Define H_A as the Hilbert completion under the positive "
                "Volterra/Sturm form A, N=ker R_global, and "
                "H_M=N cap L_M^{perp_A}."
            ),
            (
                "Let V_N be the polynomial Galerkin space and "
                "H_{M,N}=V_N cap ker R_N cap L_{M,N}^{perp_A}.  Mosco "
                "convergence H_{M,N}->H_M is the required exhaustion theorem."
            ),
            (
                "The source row is the direct Lagrange/adjoint pair "
                "E_u f=(B_P[h_u,f](s0), (P^*h_u)(s0)f(s0)).  The "
                "source-range Hardy/Green theorem constructs its A-Riesz "
                "representers uniformly, so S=E^*E is compact in the "
                "A-Hilbert source model."
            ),
            (
                "A-graph Galerkin convergence plus compact direct source-row "
                "factorization gives ||P_NSP_N-S||_{A->A}->0."
            ),
            (
                "Riesz projection convergence then gives P_{2,N}->P_2, and "
                "the normalized finite min-max tail inequality passes to "
                "H_M cap ker R_global by recovery sequences and lower limits."
            ),
        ],
        "remainingAnalyticGaps": [
            gap
            for gap in [
                None if compact_continuum_closed else compact_reason or elliptic["reason"],
                None if compact_conditional_closed else mosco_limsup["reason"],
                None if compact_conditional_closed else mosco_liminf["reason"],
                None if compact_conditional_closed else source_norm["reason"],
            ]
            if gap
        ],
        "interpretation": (
            (
                "The functional-analytic compact-source passage is proved, and "
                "the continuum trace-frame lower-bound theorem supplies the "
                "missing bounded correction right inverse.  The normalized tail "
                "estimate now passes to the full closed high block."
            )
            if tail_passes
            else (
                "The functional-analytic compact-source passage is now proved "
                "as a precise theorem.  The normalized tail estimate is not yet "
                "marked as a full continuum theorem because the finite trace-frame "
                "floor has not been promoted to a continuum uniform lower bound."
            )
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block exhaustion theorem ledger")
    print(f"  normalized epsilon_delta: {minmax['normalizedEpsilonDelta']:.12e}")
    print(f"  finite low/mid budget: {minmax['finiteLowMidSchurBudget']:.12e}")
    print(f"  abstract compact-source passage: {compact_conditional_closed}")
    print(f"  conditional min-max passage: {conditional_passage_closed}")
    print(f"  tail estimate passes to continuum: {tail_passes}")
    print(
        "  remaining: none"
        if tail_passes
        else "  remaining: continuum trace-frame lower bound"
    )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
