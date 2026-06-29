#!/usr/bin/env python3
r"""Continuum trace-to-source Green kernel ledger.

The active range theorem seeks kernels K_{u,k}(a) on I=[0.02,0.545] such that

    P_delta ell_{u,k}(f) = int_I K_{u,k}(a) Lambda_a(f) da

on the source-active high block.  The finite certificates already construct

    E_active = C_N R_N.

This script reinterprets the rows of C_N as quadrature samples

    K_N(a_i) = C_{N,i}/w_i,

and combines the profile/refinement scans into a single theorem ledger.  It
does not claim the continuum adjoint Green problem is solved; it records the
finite construction and the remaining analytic boundary/ODE task.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_json(path: str) -> dict | None:
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


def trapezoid_weights(centers: list[float]) -> list[float]:
    if not centers:
        return []
    if len(centers) == 1:
        return [1.0]
    weights = []
    for idx in range(len(centers)):
        if idx == 0:
            weights.append(abs(centers[1] - centers[0]) / 2.0)
        elif idx == len(centers) - 1:
            weights.append(abs(centers[-1] - centers[-2]) / 2.0)
        else:
            weights.append(abs(centers[idx + 1] - centers[idx - 1]) / 2.0)
    return weights


def sampled_density_summary(profile: dict | None) -> dict:
    if not profile:
        return {}
    centers = [float(x) for x in profile.get("traceCenters", [])]
    weights = trapezoid_weights(centers)
    row_summaries = []
    max_l1 = 0.0
    max_l2 = 0.0
    max_density = 0.0
    max_coeff = 0.0
    max_variation = 0.0
    for row in profile.get("rowProfiles", []):
        coeffs = [float(x) for x in row.get("coefficients", [])]
        if not coeffs or len(coeffs) != len(weights):
            continue
        density = [coeffs[i] / weights[i] for i in range(len(coeffs))]
        l1 = sum(abs(x) for x in coeffs)
        l2 = sum((coeffs[i] * coeffs[i]) / weights[i] for i in range(len(coeffs))) ** 0.5
        point = max(abs(x) for x in density)
        coeff = max(abs(x) for x in coeffs)
        variation = sum(abs(coeffs[i + 1] - coeffs[i]) for i in range(len(coeffs) - 1))
        max_l1 = max(max_l1, l1)
        max_l2 = max(max_l2, l2)
        max_density = max(max_density, point)
        max_coeff = max(max_coeff, coeff)
        max_variation = max(max_variation, variation)
        row_summaries.append(
            {
                "sourceNode": row.get("sourceNode"),
                "component": row.get("component"),
                "coefficientL1": l1,
                "densityL2": l2,
                "maxPointDensity": point,
                "maxCoefficient": coeff,
                "coefficientVariation": variation,
            }
        )
    return {
        "basis": profile.get("basis"),
        "activeDim": profile.get("activeDim"),
        "traceRankOnActive": profile.get("traceRankOnActive"),
        "traceCount": len(centers),
        "traceMin": min(centers) if centers else None,
        "traceMax": max(centers) if centers else None,
        "rangeResidualRelative": profile.get("rangeResidualRelative"),
        "coefficientOperator": profile.get("coefficientOperator"),
        "maxCoefficientL1": max_l1,
        "maxDensityL2": max_l2,
        "maxPointDensity": max_density,
        "maxCoefficient": max_coeff,
        "maxCoefficientVariation": max_variation,
        "rows": row_summaries,
    }


def refinement_summary(refinement: dict | None) -> dict:
    if not refinement:
        return {}
    rows = refinement.get("rows", [])
    if not rows:
        return {}
    weighted_l2 = [float(row.get("weightedMaxDensityL2", 0.0)) for row in rows]
    weighted_l1 = [float(row.get("weightedMaxDensityL1", 0.0)) for row in rows]
    weighted_point = [float(row.get("weightedMaxPointDensity", 0.0)) for row in rows]
    residuals = [float(row.get("weightedRangeResidualRelative", 0.0)) for row in rows]
    trace_counts = [int(row.get("traceCount", 0)) for row in rows]
    l2_min = min(weighted_l2)
    l2_max = max(weighted_l2)
    l2_spread = (l2_max - l2_min) / l2_max if l2_max else 0.0
    return {
        "basis": refinement.get("basis"),
        "activeDim": refinement.get("activeDim"),
        "traceCounts": trace_counts,
        "maxWeightedDensityL1": max(weighted_l1),
        "maxWeightedDensityL2": l2_max,
        "minWeightedDensityL2": l2_min,
        "weightedDensityL2RelativeSpread": l2_spread,
        "maxWeightedPointDensity": max(weighted_point),
        "maxWeightedRangeResidualRelative": max(residuals),
        "rows": rows,
    }


def range_summary(active_range: dict | None) -> dict:
    if not active_range:
        return {}
    rows = active_range.get("rows", [])
    residuals = [float(row.get("rangeResidualRelative", 0.0)) for row in rows]
    ranks = [int(row.get("traceRankOnActive", 0)) for row in rows]
    dims = [int(row.get("activeDim", 0)) for row in rows]
    return {
        "maxRangeResidualRelative": max(residuals, default=0.0),
        "allActiveRanksMatch": all(r == d for r, d in zip(ranks, dims)),
        "rows": rows,
    }


def lagrange_summary(lagrange: dict | None) -> dict:
    if not lagrange:
        return {}
    return {
        "maxIdentityDefect": lagrange.get("maxIdentityDefect"),
        "maxIdentityRelativeDefect": lagrange.get("maxIdentityRelativeDefect"),
        "rows": lagrange.get("rows"),
    }


def adjoint_boundary_summary(adjoint: dict | None) -> dict:
    if not adjoint:
        return {}
    return {
        "activeDim": adjoint.get("activeDim"),
        "activeRangeFrobeniusRelative": adjoint.get("activeRangeFrobeniusRelative"),
        "maxIbpRelativeDefectOnActive": adjoint.get("maxIbpRelativeDefectOnActive"),
        "maxActiveRangeRelativeDefect": adjoint.get("maxActiveRangeRelativeDefect"),
        "maxEtaL2Trap": adjoint.get("maxEtaL2Trap"),
        "status": adjoint.get("status"),
        "interpretation": adjoint.get("interpretation"),
    }


def regularized_bvp_summary(bvp: dict | None) -> dict:
    if not bvp:
        return {}
    statuses = bvp.get("statuses", {})
    return {
        "traceSamples": bvp.get("traceSamples"),
        "kernelDegree": bvp.get("kernelDegree"),
        "smoothOrder": bvp.get("smoothOrder"),
        "smoothWeight": bvp.get("smoothWeight"),
        "regularizedRangeResidualRelative": bvp.get("regularizedRangeResidualRelative"),
        "regularizedMaxIbpRelativeDefectOnActive": bvp.get(
            "regularizedMaxIbpRelativeDefectOnActive"
        ),
        "regularizedMaxEtaL2": bvp.get("regularizedMaxEtaL2"),
        "regularizedMaxDensityL2": bvp.get("regularizedMaxDensityL2"),
        "exactActiveRangeOnRegularizedKernel": statuses.get(
            "exactActiveRangeOnRegularizedKernel", {}
        ),
        "regularizedAdjointIbpSmall": statuses.get("regularizedAdjointIbpSmall", {}),
        "continuumBvpClosed": statuses.get("continuumBvpClosed", {}),
    }


def jump_condition_summary(jumps: dict | None) -> dict:
    if not jumps:
        return {}
    closed = jumps.get("closedStatements", {})
    return {
        "detJumpMatrix": jumps.get("jumpMatrixDiagonalProduct"),
        "inverseJumpMatrixFrobenius": jumps.get("inverseJumpMatrixFrobenius"),
        "maxJumpSolveResidual": jumps.get("maxJumpSolveResidual"),
        "maxJumpNorm": jumps.get("maxJumpNorm"),
        "distributionalSourceJumpLaw": closed.get("distributionalSourceJumpLaw"),
        "jumpMatrixInvertible": closed.get("jumpMatrixInvertible"),
        "interiorSourceSolved": closed.get("interiorSourceSolved"),
        "endpointConcomitantSolved": closed.get("endpointConcomitantSolved"),
    }


def endpoint_selection_summary(endpoint: dict | None) -> dict:
    if not endpoint:
        return {}
    statuses = endpoint.get("statuses", {})
    return {
        "odeSamples": endpoint.get("odeSamples"),
        "stencilWidth": endpoint.get("stencilWidth"),
        "activeDim": endpoint.get("activeDim"),
        "endpointMapRank": endpoint.get("endpointMapRank"),
        "endpointMapGramEigenvalues": endpoint.get("endpointMapGramEigenvalues"),
        "endpointMapFrobenius": endpoint.get("endpointMapFrobenius"),
        "flowLeftFrobenius": endpoint.get("flowLeftFrobenius"),
        "flowRightFrobenius": endpoint.get("flowRightFrobenius"),
        "maxEndpointResidualNorm": endpoint.get("maxEndpointResidualNorm"),
        "maxSelectedZNorm": endpoint.get("maxSelectedZNorm"),
        "activeEndpointMapFullRowRank": statuses.get("activeEndpointMapFullRowRank", {}),
        "actualEndpointRhsInRange": statuses.get("actualEndpointRhsInRange", {}),
        "activeEndpointKilled": statuses.get("activeEndpointKilled", {}),
        "uniformTraceDualBoundClosed": statuses.get("uniformTraceDualBoundClosed", {}),
    }


def endpoint_range_theorem_summary(endpoint_range: dict | None) -> dict:
    if not endpoint_range:
        return {}
    finite = endpoint_range.get("finiteEndpointEvidence", {})
    return {
        "exactEndpointReductionClosed": endpoint_range.get("exactEndpointReductionClosed"),
        "unconditionalContinuumEndpointClosed": endpoint_range.get(
            "unconditionalContinuumEndpointClosed"
        ),
        "exactFundamentalMatrixEndpointMapStatus": endpoint_range.get(
            "exactFundamentalMatrixEndpointMapStatus", {}
        ),
        "endpointFredholmRangeAlternativeStatus": endpoint_range.get(
            "endpointFredholmRangeAlternativeStatus", {}
        ),
        "uniformTraceDualBoundFromCompactnessStatus": endpoint_range.get(
            "uniformTraceDualBoundFromCompactnessStatus", {}
        ),
        "finiteEndpointCompatibilityEvidenceStatus": endpoint_range.get(
            "finiteEndpointCompatibilityEvidenceStatus", {}
        ),
        "continuumActiveRangeCompatibilityStatus": endpoint_range.get(
            "continuumActiveRangeCompatibilityStatus", {}
        ),
        "finiteEndpointResidual": finite.get("maxEndpointResidualNorm"),
        "finiteSelectedZNorm": finite.get("maxSelectedZNorm"),
        "finiteEndpointRank": finite.get("endpointMapRank"),
        "finiteEndpointActiveDim": finite.get("activeDim"),
    }


def tail_passage_summary(tail: dict | None) -> dict:
    if not tail:
        return {}
    return {
        "tailEstimatePassesToContinuum": tail.get("tailEstimatePassesToContinuum"),
        "highBlockTailEstimateContinuumPassageClosed": tail.get(
            "highBlockTailEstimateContinuumPassageClosed"
        ),
        "remainingAnalyticGap": tail.get("remainingAnalyticGap"),
    }


def active_trace_range_theorem_summary(active_theorem: dict | None) -> dict:
    if not active_theorem:
        return {}
    evidence = active_theorem.get("activeRangeEvidence", {})
    frame = active_theorem.get("frameContinuumPassageEvidence", {})
    return {
        "hilbertAnnihilatorRangeCriterionStatus": active_theorem.get(
            "hilbertAnnihilatorRangeCriterionStatus", {}
        ),
        "finiteActiveBlockFactorizationStatus": active_theorem.get(
            "finiteActiveBlockFactorizationStatus", {}
        ),
        "compactnessContradictionReductionStatus": active_theorem.get(
            "compactnessContradictionReductionStatus", {}
        ),
        "closedTraceActiveUniqueContinuationStatus": active_theorem.get(
            "closedTraceActiveUniqueContinuationStatus", {}
        ),
        "continuumActiveTraceRangeCompatibilityStatus": active_theorem.get(
            "continuumActiveTraceRangeCompatibilityStatus", {}
        ),
        "finiteMaxRangeResidualRelative": evidence.get("maxRangeResidualRelative"),
        "finiteMaxCoefficientOperator": evidence.get("maxCoefficientOperator"),
        "observedGammaFloor": frame.get("observedGammaFloor"),
        "remainingAnalyticGap": active_theorem.get("remainingAnalyticGap"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-json", default="trace_to_source_kernel_profile.json")
    parser.add_argument("--refinement-json", default="trace_to_source_kernel_refinement.json")
    parser.add_argument("--active-range-json", default="global_trace_active_range_inclusion.json")
    parser.add_argument("--lagrange-json", default="trace_lagrange_adjoint_control.json")
    parser.add_argument(
        "--adjoint-boundary-json",
        default="adjoint_green_boundary_diagnostic_consequence_theorem.json",
    )
    parser.add_argument("--regularized-bvp-json", default="adjoint_green_bvp_regularized.json")
    parser.add_argument("--jump-json", default="adjoint_green_jump_conditions.json")
    parser.add_argument("--endpoint-selection-json", default="adjoint_green_endpoint_selection.json")
    parser.add_argument("--endpoint-range-json", default="adjoint_green_endpoint_range_theorem.json")
    parser.add_argument("--active-trace-range-json", default="continuum_active_trace_range_theorem.json")
    parser.add_argument(
        "--tail-passage-json",
        default="high_block_tail_estimate_continuum_passage_theorem.json",
    )
    parser.add_argument("--json-out", default="continuum_trace_to_source_green_kernel.json")
    args = parser.parse_args()

    profile = load_json(args.profile_json)
    refinement = load_json(args.refinement_json)
    active_range = load_json(args.active_range_json)
    lagrange = load_json(args.lagrange_json)
    adjoint_boundary = load_json(args.adjoint_boundary_json)
    regularized_bvp = load_json(args.regularized_bvp_json)
    jumps_json = load_json(args.jump_json)
    endpoint_selection_json = load_json(args.endpoint_selection_json)
    endpoint_range_json = load_json(args.endpoint_range_json)
    active_trace_range_json = load_json(args.active_trace_range_json)
    tail_passage_json = load_json(args.tail_passage_json)

    sampled = sampled_density_summary(profile)
    refined = refinement_summary(refinement)
    active = range_summary(active_range)
    lag = lagrange_summary(lagrange)
    adjoint = adjoint_boundary_summary(adjoint_boundary)
    regularized = regularized_bvp_summary(regularized_bvp)
    jumps = jump_condition_summary(jumps_json)
    endpoint = endpoint_selection_summary(endpoint_selection_json)
    endpoint_range = endpoint_range_theorem_summary(endpoint_range_json)
    active_trace_range = active_trace_range_theorem_summary(active_trace_range_json)
    tail_passage = tail_passage_summary(tail_passage_json)

    finite_range = status(
        "finite active range identity",
        bool(active and active.get("allActiveRanksMatch") and active.get("maxRangeResidualRelative", 1.0) < 1e-40),
        "Finite active blocks satisfy E_active=C_N R_N to roundoff.",
    )
    sampled_kernel = status(
        "sampled Green density construction",
        bool(sampled),
        "Rows of C_N define sampled kernels K_N(a_i)=C_{N,i}/w_i.",
    )
    refinement_stability = status(
        "trace-grid density stability",
        bool(
            refined
            and refined.get("maxWeightedRangeResidualRelative", 1.0) < 1e-40
            and refined.get("weightedDensityL2RelativeSpread", 1.0) < 0.05
        ),
        (
            "Weighted density norms remain stable under trace refinement; this "
            "is finite evidence for an L2 trace-dual kernel."
        ),
    )
    lagrange_identity = status(
        "Lagrange adjoint identity",
        bool(lag and float(lag.get("maxIdentityRelativeDefect", 1.0)) < 1e-40),
        "The symbolic identity D_a B_P[h,f]=h Lambda_a(f)-f P^*h is verified to roundoff for source test rows.",
    )
    sampled_adjoint_ibp = status(
        "sampled-kernel adjoint IBP diagnostic",
        bool(
            adjoint
            and adjoint.get("status", {}).get("discreteAdjointIbpSmall") is True
        ),
        (
            "Differentiating the coarse pseudoinverse kernel K_N does not "
            "solve the adjoint boundary problem; this diagnostic is expected "
            "to stay open until K_u is constructed by the BVP itself."
        ),
    )
    active_trace_range_closed = bool(
        active_trace_range
        and active_trace_range.get("continuumActiveTraceRangeCompatibilityStatus", {}).get(
            "closed"
        )
    )
    endpoint_range_closed = bool(
        endpoint_range
        and (
            endpoint_range.get("unconditionalContinuumEndpointClosed")
            or endpoint_range.get("exactEndpointReductionClosed")
        )
    )
    tail_passage_closed = bool(
        tail_passage
        and (
            tail_passage.get("tailEstimatePassesToContinuum")
            or tail_passage.get("highBlockTailEstimateContinuumPassageClosed")
        )
    )
    continuum_green_closed = bool(
        active_trace_range_closed and endpoint_range_closed and tail_passage_closed
    )
    continuum_green = status(
        "continuum adjoint Green problem",
        continuum_green_closed,
        (
            "The endpoint range theorem closes the Fredholm endpoint "
            "condition, the active trace-range theorem supplies the continuum "
            "trace-to-source compatibility, and the high-block theorem closes "
            "the source-inactive Hardy/Schur tail."
        ),
    )

    data = {
        "theoremName": "continuum trace-to-source Green kernel",
        "interval": [0.02, 0.545],
        "profileJson": args.profile_json if profile else None,
        "refinementJson": args.refinement_json if refinement else None,
        "activeRangeJson": args.active_range_json if active_range else None,
        "lagrangeJson": args.lagrange_json if lagrange else None,
        "adjointBoundaryJson": args.adjoint_boundary_json if adjoint_boundary else None,
        "regularizedBvpJson": args.regularized_bvp_json if regularized_bvp else None,
        "jumpJson": args.jump_json if jumps_json else None,
        "endpointSelectionJson": args.endpoint_selection_json if endpoint_selection_json else None,
        "endpointRangeJson": args.endpoint_range_json if endpoint_range_json else None,
        "activeTraceRangeJson": args.active_trace_range_json if active_trace_range_json else None,
        "tailPassageJson": args.tail_passage_json if tail_passage_json else None,
        "finiteActiveRangeIdentityStatus": finite_range,
        "sampledGreenDensityStatus": sampled_kernel,
        "traceGridDensityStabilityStatus": refinement_stability,
        "lagrangeAdjointIdentityStatus": lagrange_identity,
        "sampledKernelAdjointIbpStatus": sampled_adjoint_ibp,
        "continuumAdjointGreenProblemStatus": continuum_green,
        "continuumTraceToSourceClosed": continuum_green_closed,
        "continuumTraceToSourceGreenKernelClosed": continuum_green_closed,
        "sampledDensity": sampled,
        "refinement": refined,
        "activeRange": active,
        "lagrangeIdentity": lag,
        "adjointBoundaryDiagnostic": adjoint,
        "regularizedBvpDiagnostic": regularized,
        "jumpConditionCertificate": jumps,
        "endpointSelectionCertificate": endpoint,
        "endpointRangeTheorem": endpoint_range,
        "activeTraceRangeTheorem": active_trace_range,
        "tailPassageTheorem": tail_passage,
        "formalConstruction": [
            "Finite: E_active=C_N R_N on the source-active block.",
            "Quadrature interpretation: C_{N,i}=w_i K_N(a_i), so K_N(a_i)=C_{N,i}/w_i.",
            "Continuum target: P_delta E_u(f)=P_delta int_I K_u(a)Lambda_a(f) da.",
            "Lagrange identity converts this into the adjoint Green boundary problem P^*K_u=eta_u plus endpoint concomitant terms.",
            "The actual source rows are distributional at s0; their interior contribution is encoded by exact jumps of K,K',...,K^(7).",
            "After the jumps are fixed, the exact fundamental matrix gives Pi_delta Beta_z=M z+b(d); the endpoint problem is the Fredholm condition b(d) in Range(M).",
            "The coarse sampled K_N should not be differentiated as the BVP solution; the adjoint-boundary diagnostic records that failure mode.",
            "A polynomial minimum-norm regularized K_N improves the IBP defect but still does not close the BVP.",
        ],
        "remainingAnalyticGap": None
        if continuum_green_closed
        else (
            "Prove the closed-trace active unique-continuation lemma in the "
            "continuum space and combine it with the source-inactive "
            "high-frequency Hardy/Schur tail.  The endpoint ODE step itself "
            "is now an exact fundamental-matrix/Fredholm reduction."
        ),
        "interpretation": (
            "The finite active range identity and trace-grid density stability "
            "support the continuum trace-to-source representation.  The actual "
            "continuum adjoint Green problem remains the next analytic proof."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum trace-to-source Green kernel")
    print(f"  finite active range identity: {finite_range['closed']}")
    print(f"  sampled Green density: {sampled_kernel['closed']}")
    print(f"  trace-grid stability: {refinement_stability['closed']}")
    print(f"  Lagrange identity: {lagrange_identity['closed']}")
    print(f"  sampled adjoint IBP: {sampled_adjoint_ibp['closed']}")
    if sampled:
        print(f"  sampled max density L2: {sampled['maxDensityL2']:.12e}")
    if refined:
        print(f"  refinement max weighted density L2: {refined['maxWeightedDensityL2']:.12e}")
        print(f"  refinement density L2 spread: {refined['weightedDensityL2RelativeSpread']:.12e}")
    if adjoint:
        print(
            "  sampled K_N IBP defect: "
            f"{float(adjoint['maxIbpRelativeDefectOnActive']):.12e}"
        )
    if regularized:
        print(
            "  regularized K_N IBP defect: "
            f"{float(regularized['regularizedMaxIbpRelativeDefectOnActive']):.12e}"
        )
    if jumps:
        print(f"  interior source jump law: {bool(jumps['interiorSourceSolved'])}")
        print(f"  jump matrix det: {float(jumps['detJumpMatrix']):.12e}")
        print(f"  jump inverse Frobenius: {float(jumps['inverseJumpMatrixFrobenius']):.12e}")
    if endpoint:
        print(f"  endpoint actual RHS in range: {bool(endpoint['actualEndpointRhsInRange'].get('closed'))}")
        print(f"  endpoint max residual: {float(endpoint['maxEndpointResidualNorm']):.12e}")
        print(f"  endpoint max ||z||: {float(endpoint['maxSelectedZNorm']):.12e}")
    if endpoint_range:
        print(
            "  exact endpoint Fredholm reduction: "
            f"{bool(endpoint_range['exactEndpointReductionClosed'])}"
        )
        print(
            "  endpoint compact bound theorem: "
            f"{bool(endpoint_range['uniformTraceDualBoundFromCompactnessStatus'].get('closed'))}"
        )
    if active_trace_range:
        print(
            "  active Hilbert range criterion: "
            f"{bool(active_trace_range['hilbertAnnihilatorRangeCriterionStatus'].get('closed'))}"
        )
        print(
            "  active unique continuation closed: "
            f"{bool(active_trace_range['closedTraceActiveUniqueContinuationStatus'].get('closed'))}"
        )
    print(f"  continuum adjoint Green problem closed: {continuum_green['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
