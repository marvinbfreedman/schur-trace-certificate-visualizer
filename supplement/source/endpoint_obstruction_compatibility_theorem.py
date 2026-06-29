#!/usr/bin/env python3
r"""Endpoint obstruction compatibility theorem ledger.

For the exact endpoint equation

    M z = -b(d),

compatibility is

    w^* b(d)=0 for every w in ker M^*.

This ledger records the non-circular Fredholm obstruction formula.  A left
endpoint obstruction is not merely a sampled vector; in the exact boundary
problem it corresponds, through Green's formula, to a primal homogeneous
endpoint mode.  For any jump-source Green coefficient K_d,

    w^* b(d) = int_I K_d(a) P f_w(a) da - d(f_w),

where f_w is the primal endpoint test represented by the obstruction.  Hence
endpoint compatibility follows if the Volterra/Sturm primal obstruction
theorem holds:

    P f_w = 0 and endpoint-dual boundary conditions  =>  d(f_w)=0

for every actual active source row d=d_u.

There is a sharper alternative: if the exact endpoint map M has full active row
rank, then ker(M^*)={0} and compatibility is vacuous.  The companion rank-ball
certificate records the quantitative singular margin needed to close that
route.
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


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def endpoint_evidence(endpoint: dict | None) -> dict:
    if not endpoint:
        return {}
    statuses = endpoint.get("statuses", {})
    return {
        "activeDim": endpoint.get("activeDim"),
        "endpointMapRank": endpoint.get("endpointMapRank"),
        "endpointMapGramEigenvalues": endpoint.get("endpointMapGramEigenvalues"),
        "maxEndpointResidualNorm": endpoint.get("maxEndpointResidualNorm"),
        "maxSelectedZNorm": endpoint.get("maxSelectedZNorm"),
        "leftObstructionVector": endpoint.get("leftObstructionVector"),
        "leftObstructionResidualNorm": endpoint.get("leftObstructionResidualNorm"),
        "maxLeftObstructionPairing": endpoint.get("maxLeftObstructionPairing"),
        "maxLeftObstructionRelativePairing": endpoint.get(
            "maxLeftObstructionRelativePairing"
        ),
        "actualEndpointRhsInRange": statuses.get("actualEndpointRhsInRange", {}),
        "activeEndpointKilled": statuses.get("activeEndpointKilled", {}),
        "finiteLeftObstructionAnnihilatesActualRows": statuses.get(
            "finiteLeftObstructionAnnihilatesActualRows", {}
        ),
        "interpretation": (
            "Finite sampled-flow endpoint evidence.  A tiny full solve "
            "residual only proves the sampled endpoint equation can be solved; "
            "it does not by itself prove annihilation of continuum Fredholm "
            "obstructions."
        ),
    }


def rank_evidence(rank: dict | None) -> dict:
    if not rank:
        return {}
    return {
        "conditionalBallFullRowRankStatus": rank.get(
            "conditionalBallFullRowRankStatus", {}
        ),
        "rigorousEndpointMapFullRowRankStatus": rank.get(
            "rigorousEndpointMapFullRowRankStatus", {}
        ),
        "kernelMStarZeroStatus": rank.get("kernelMStarZeroStatus", {}),
        "centerSigmaMin": rank.get("centerSigmaMin"),
        "centerSigmaMax": rank.get("centerSigmaMax"),
        "lowerSigmaMin": rank.get("lowerSigmaMin"),
        "pinvNormUpperBound": rank.get("pinvNormUpperBound"),
        "maxAllowedFrobeniusRadiusForRank": rank.get(
            "maxAllowedFrobeniusRadiusForRank"
        ),
        "maxAllowedRelativeFrobeniusRadiusForRank": rank.get(
            "maxAllowedRelativeFrobeniusRadiusForRank"
        ),
        "remainingAnalyticGap": rank.get("remainingAnalyticGap"),
    }


def center_stability_evidence(center: dict | None) -> dict:
    if not center:
        return {}
    return {
        "sampledCenterStabilityStatus": center.get(
            "sampledCenterStabilityStatus", {}
        ),
        "exactSegmentCenterStabilityStatus": center.get(
            "exactSegmentCenterStabilityStatus", {}
        ),
        "rankBallCenterReadyStatus": center.get("rankBallCenterReadyStatus", {}),
        "sampledFrobeniusSpread": center.get("sampledFrobeniusSpread"),
        "exactSegmentFrobeniusSpread": center.get("exactSegmentFrobeniusSpread"),
        "conclusion": center.get("conclusion"),
        "nextTarget": center.get("nextTarget"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-selection-json", default="adjoint_green_endpoint_selection.json")
    parser.add_argument("--endpoint-range-json", default="adjoint_green_endpoint_range_theorem.json")
    parser.add_argument("--jump-json", default="adjoint_green_jump_conditions.json")
    parser.add_argument("--lagrange-json", default="trace_lagrange_adjoint_control.json")
    parser.add_argument("--rank-json", default="endpoint_map_rank_ball_certificate.json")
    parser.add_argument("--center-stability-json", default="endpoint_map_center_stability.json")
    parser.add_argument("--json-out", default="endpoint_obstruction_compatibility_theorem.json")
    args = parser.parse_args()

    endpoint = load_optional(args.endpoint_selection_json)
    endpoint_range = load_optional(args.endpoint_range_json)
    jumps = load_optional(args.jump_json)
    lagrange = load_optional(args.lagrange_json)
    rank = load_optional(args.rank_json)
    center = load_optional(args.center_stability_json)

    finite = endpoint_evidence(endpoint)
    rank_summary = rank_evidence(rank)
    center_summary = center_stability_evidence(center)
    lagrange_closed = bool(
        lagrange and float(lagrange.get("maxIdentityRelativeDefect", 1.0)) < 1e-40
    )
    jumps_closed = bool(
        jumps
        and jumps.get("closedStatements", {}).get("interiorSourceSolved")
        and jumps.get("closedStatements", {}).get("jumpMatrixInvertible")
    )
    endpoint_reduction_closed = bool(
        endpoint_range
        and endpoint_range.get("exactFundamentalMatrixEndpointMapStatus", {}).get("closed")
        and endpoint_range.get("endpointFredholmRangeAlternativeStatus", {}).get("closed")
    )
    if finite and finite.get("finiteLeftObstructionAnnihilatesActualRows"):
        finite_annihilation = bool(
            finite["finiteLeftObstructionAnnihilatesActualRows"].get("closed")
        )
    else:
        finite_annihilation = False
    rigorous_full_rank = bool(
        rank_summary.get("rigorousEndpointMapFullRowRankStatus", {}).get("closed")
        or (
            endpoint_range
            and endpoint_range.get("endpointFullActiveRowRankStatus", {}).get(
                "closed"
            )
        )
        or (
            endpoint_range
            and endpoint_range.get("unconditionalContinuumEndpointClosed")
        )
    )

    fredholm_obstruction = status(
        "endpoint Fredholm obstruction formula",
        endpoint_reduction_closed,
        (
            "The exact endpoint problem is Mz=-b(d); compatibility is "
            "w^*b(d)=0 for every w in ker M^*."
        ),
    )
    green_pairing = status(
        "Green pairing formula for obstruction",
        lagrange_closed and jumps_closed and endpoint_reduction_closed,
        (
            "Green's identity and the exact jump law give "
            "w^*b(d)=int K_d P f_w - d(f_w) for the primal endpoint test f_w."
        ),
    )
    finite_status = status(
        "finite sampled near-obstruction annihilation",
        finite_annihilation,
        (
            "The resolved relative-threshold left endpoint direction must "
            "annihilate the actual source endpoint vectors.  This is stronger "
            "than solving the finite endpoint equation with the full sampled "
            "matrix."
        ),
    )
    primal_obstruction = status(
        "Volterra/Sturm primal obstruction annihilation",
        False,
        (
            "Still need the non-circular continuum theorem: every primal "
            "homogeneous endpoint obstruction f_w annihilates the actual "
            "active source rows d_u."
        ),
    )
    full_rank_vacuity = status(
        "full endpoint row rank vacuity",
        rigorous_full_rank,
        (
            "If the endpoint-range interval theorem proves full active row "
            "rank of M, then ker(M^*)={0}; there is no Fredholm obstruction "
            "to annihilate.  Older rank-ball files are retained only as "
            "diagnostic evidence."
        ),
    )
    theorem_closed = bool(
        green_pairing["closed"]
        and (full_rank_vacuity["closed"] or primal_obstruction["closed"])
    )
    compatibility = status(
        "continuum endpoint obstruction compatibility",
        theorem_closed,
        (
            "Closed either by rigorous full active row rank of M, which makes "
            "ker(M^*) trivial, or by combining the Green pairing formula with "
            "the Volterra/Sturm primal obstruction annihilation theorem."
        ),
    )

    data = {
        "theoremName": "endpoint obstruction compatibility",
        "endpointSelectionJson": args.endpoint_selection_json if endpoint else None,
        "endpointRangeJson": args.endpoint_range_json if endpoint_range else None,
        "jumpJson": args.jump_json if jumps else None,
        "lagrangeJson": args.lagrange_json if lagrange else None,
        "rankJson": args.rank_json if rank else None,
        "centerStabilityJson": args.center_stability_json if center else None,
        "endpointFredholmObstructionFormulaStatus": fredholm_obstruction,
        "greenPairingFormulaForObstructionStatus": green_pairing,
        "finiteSampledObstructionAnnihilationStatus": finite_status,
        "fullEndpointRowRankVacuityStatus": full_rank_vacuity,
        "volterraSturmPrimalObstructionAnnihilationStatus": primal_obstruction,
        "continuumEndpointObstructionCompatibilityStatus": compatibility,
        "continuumEndpointObstructionCompatibilityClosed": theorem_closed,
        "finiteEvidence": finite,
        "rankBallEvidence": rank_summary,
        "endpointRangeFullRankEvidence": (
            endpoint_range.get("endpointFullActiveRowRankStatus", {})
            if endpoint_range
            else {}
        ),
        "centerStabilityEvidence": center_summary,
        "proof": [
            (
                "For the exact fundamental matrix endpoint map, the endpoint "
                "equation is Mz=-b(d).  The finite-dimensional Fredholm "
                "alternative gives compatibility iff ker(M^*) annihilates b(d)."
            ),
            (
                "Let w in ker(M^*) and let f_w be the corresponding primal "
                "endpoint test.  The condition w^*M=0 says all homogeneous "
                "adjoint endpoint concomitants vanish against f_w."
            ),
            (
                "For the jump Green coefficient K_d, integrate "
                "dB_P[K_d,f_w]=K_d P f_w-f_w P^*K_d across I.  The jump law "
                "turns P^*K_d into the source row d."
            ),
            (
                "Thus w^*b(d)=int_I K_d(a)P f_w(a)da-d(f_w).  If f_w is a "
                "primal homogeneous endpoint obstruction, this reduces to "
                "w^*b(d)=-d(f_w)."
            ),
            (
                "The remaining continuum theorem is therefore d_u(f_w)=0 for "
                "the actual active source family and every primal homogeneous "
                "endpoint obstruction f_w."
            ),
            (
                "Alternatively, a rigorous rank-ball proof that sigma_min(M)>0 "
                "makes ker(M^*) trivial and closes endpoint compatibility "
                "without an obstruction-annihilation theorem."
            ),
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else (
            center_summary.get("nextTarget")
            or rank_summary.get("remainingAnalyticGap")
            or (
                "Either prove the Volterra/Sturm primal obstruction "
                "annihilation theorem, or prove rigorous full active row rank "
                "of the exact endpoint map M."
            )
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Endpoint obstruction compatibility")
    print(f"  Fredholm obstruction formula: {fredholm_obstruction['closed']}")
    print(f"  Green pairing formula: {green_pairing['closed']}")
    print(f"  finite sampled annihilation: {finite_status['closed']}")
    print(f"  full row-rank vacuity: {full_rank_vacuity['closed']}")
    if finite:
        print(
            "  finite max endpoint residual: "
            f"{float(finite['maxEndpointResidualNorm']):.12e}"
        )
        if finite.get("maxLeftObstructionRelativePairing") is not None:
            print(
                "  finite max relative left-obstruction pairing: "
                f"{float(finite['maxLeftObstructionRelativePairing']):.12e}"
            )
    print(f"  primal obstruction annihilation closed: {primal_obstruction['closed']}")
    print(f"  continuum compatibility closed: {compatibility['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
