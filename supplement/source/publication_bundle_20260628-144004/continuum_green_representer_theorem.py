#!/usr/bin/env python3
r"""Continuum Green representer construction ledger.

The valid source-range Hardy/Green target is the direct source-row theorem.
On the completed high closed-trace A-space H_hi, define

    E_u f = (ell_{u,1}(f), ell_{u,2}(f))
          = (B_P[h_u,f](s0), (P^*h_u)(s0) f(s0)).

The desired Green representers are the A-Riesz representatives

    ell_{u,k}(f) = <g_{u,k}, f>_A,        k=1,2.

Their 2x2 A-Gram is uniformly bounded on u in [0.08,0.52] exactly when the
source rows are uniformly A-bounded.  The older fixed-jet construction

    g_{u,1}=sum_j b_j(u)k_j^hi,   g_{u,2}=p(u)k_0^hi

is only a sufficient route.  Existing scans show those raw fixed-jet constants
grow, so the proof should instead use active trace range inclusion plus the
source-inactive Hardy/Schur tail estimate for the actual source rows.
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


def fixed_summary(fixed: dict | None) -> dict:
    if not fixed:
        return {}
    rows = fixed.get("rows", [])
    return {
        "bases": fixed.get("bases"),
        "cutoff": fixed.get("cutoff"),
        "coefficientEnvelope": fixed.get("coefficientEnvelope"),
        "rows": rows,
        "maxEvalRepresenterNorm2": max(
            (float(row.get("evalRepresenterNorm2", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxEvalCoverNorm2": max(
            (float(row.get("evalCoverNorm2", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxBoundaryCoverNorm2": max(
            (float(row.get("boundaryCoverNorm2", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxBoundaryDerivativeCoverNorm2": max(
            (float(row.get("boundaryDerivativeCoverNorm2", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxJetRangeRelativeDefect": max(
            (float(row.get("maxJetRangeRelativeDefect", 0.0)) for row in rows),
            default=0.0,
        ),
    }


def adjoint_summary(adjoint: dict | None) -> dict:
    if not adjoint:
        return {}
    return {
        "basis": adjoint.get("basis"),
        "constraints": adjoint.get("constraints"),
        "cutoff": adjoint.get("cutoff"),
        "evalRepresenterNorm2": adjoint.get("evalRepresenterNorm2"),
        "evalRangeRelativeDefect": adjoint.get("evalRangeRelativeDefect"),
        "maxFactorRelativeDefect": adjoint.get("maxFactorRelativeDefect"),
        "coverEvalNorm2": adjoint.get("coverEvalNorm2"),
        "coverDEvalNorm2": adjoint.get("coverDEvalNorm2"),
    }


def boundary_summary(boundary: dict | None) -> dict:
    if not boundary:
        return {}
    return {
        "basis": boundary.get("basis"),
        "constraints": boundary.get("constraints"),
        "cutoff": boundary.get("cutoff"),
        "maxJetRangeRelativeDefect": boundary.get("maxJetRangeRelativeDefect"),
        "maxFactorRelativeDefect": boundary.get("maxFactorRelativeDefect"),
        "coverBoundaryNorm2": boundary.get("coverBoundaryNorm2"),
        "coverDBoundaryNorm2": boundary.get("coverDBoundaryNorm2"),
    }


def closed_green_summary(closed_green: dict | None) -> dict:
    if not closed_green:
        return {}
    rows = closed_green.get("rows", [])
    return {
        "basis": closed_green.get("basis"),
        "constraints": closed_green.get("constraints"),
        "cutoff": closed_green.get("cutoff"),
        "highModes": closed_green.get("highModes"),
        "lambdaMinHigh": closed_green.get("lambdaMinHigh"),
        "maxRangeRelativeDefect": closed_green.get("maxRangeRelativeDefect"),
        "maxFiniteEConstant": max(
            (float(row.get("EConstant", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxFiniteDEConstant": max(
            (float(row.get("dEConstant", 0.0)) for row in rows),
            default=0.0,
        ),
        "rows": rows,
    }


def global_observability_summary(global_obs: dict | None) -> dict:
    if not global_obs:
        return {}
    rows = global_obs.get("rows", [])
    return {
        "bases": global_obs.get("bases"),
        "sourceGrid": global_obs.get("sourceGrid"),
        "constraintRule": global_obs.get("constraintRule"),
        "constraintRatio": global_obs.get("constraintRatio"),
        "maxHighFullFrac": max(
            (float(row.get("maxHighFullFrac", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxDerivativeHighFullFrac": max(
            (float(row.get("maxDerivativeHighFullFrac", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxHighConstant": max(
            (float(row.get("maxHighConstant", 0.0)) for row in rows),
            default=0.0,
        ),
        "maxDerivativeHighConstant": max(
            (float(row.get("maxDerivativeHighConstant", 0.0)) for row in rows),
            default=0.0,
        ),
        "rows": rows,
    }


def trace_to_source_summary(trace_to_source: dict | None) -> dict:
    if not trace_to_source:
        return {}
    l1_values = []
    for row in trace_to_source.get("rowProfiles", []):
        row_stats = row.get("stats", {})
        l1_values.append(float(row_stats.get("l1Grid", 0.0)))
    return {
        "basis": trace_to_source.get("basis"),
        "activeDim": trace_to_source.get("activeDim"),
        "traceRankOnActive": trace_to_source.get("traceRankOnActive"),
        "coefficientOperator": trace_to_source.get("coefficientOperator"),
        "coefficientFrobenius": trace_to_source.get("coefficientFrobenius"),
        "rangeResidualRelative": trace_to_source.get("rangeResidualRelative"),
        "maxKernelL1Grid": max(l1_values, default=0.0),
    }


def inactive_tail_summary(inactive_tail: dict | None) -> dict:
    if not inactive_tail:
        return {}
    return {
        "sourceGrid": inactive_tail.get("sourceGrid"),
        "activeCount": inactive_tail.get("activeCount"),
        "continuumErrorBound": inactive_tail.get("continuumErrorBound"),
        "worstContinuumInactiveFracUpper": inactive_tail.get(
            "worstContinuumInactiveFracUpper"
        ),
        "worstActiveTraceKernelSourceFrac": inactive_tail.get(
            "worstActiveTraceKernelSourceFrac"
        ),
        "worstFullTraceKernelSourceFrac": inactive_tail.get(
            "worstFullTraceKernelSourceFrac"
        ),
        "sampledActiveKernelExclusionPasses": inactive_tail.get(
            "sampledActiveKernelExclusionPasses"
        ),
        "finiteBudgetDiagnosticPasses": inactive_tail.get("finiteBudgetDiagnosticPasses"),
        "globalSchurTailTheoremClosed": inactive_tail.get("globalSchurTailTheoremClosed"),
    }


def trace_green_summary(trace_green: dict | None) -> dict:
    if not trace_green:
        return {}
    sampled = trace_green.get("sampledDensity", {})
    refined = trace_green.get("refinement", {})
    regularized = trace_green.get("regularizedBvpDiagnostic", {})
    jumps = trace_green.get("jumpConditionCertificate", {})
    endpoint = trace_green.get("endpointSelectionCertificate", {})
    endpoint_range = trace_green.get("endpointRangeTheorem", {})
    active_range_theorem = trace_green.get("activeTraceRangeTheorem", {})
    return {
        "finiteActiveRangeIdentityStatus": trace_green.get(
            "finiteActiveRangeIdentityStatus"
        ),
        "sampledGreenDensityStatus": trace_green.get("sampledGreenDensityStatus"),
        "traceGridDensityStabilityStatus": trace_green.get(
            "traceGridDensityStabilityStatus"
        ),
        "lagrangeAdjointIdentityStatus": trace_green.get("lagrangeAdjointIdentityStatus"),
        "continuumAdjointGreenProblemStatus": trace_green.get(
            "continuumAdjointGreenProblemStatus"
        ),
        "sampledKernelAdjointIbpStatus": trace_green.get(
            "sampledKernelAdjointIbpStatus"
        ),
        "sampledMaxDensityL2": sampled.get("maxDensityL2"),
        "sampledMaxPointDensity": sampled.get("maxPointDensity"),
        "refinementMaxWeightedDensityL2": refined.get("maxWeightedDensityL2"),
        "refinementDensityL2RelativeSpread": refined.get(
            "weightedDensityL2RelativeSpread"
        ),
        "remainingAnalyticGap": trace_green.get("remainingAnalyticGap"),
        "adjointBoundaryDiagnostic": trace_green.get("adjointBoundaryDiagnostic"),
        "regularizedBvpDiagnostic": regularized,
        "jumpConditionCertificate": jumps,
        "endpointSelectionCertificate": endpoint,
        "endpointRangeTheorem": endpoint_range,
        "activeTraceRangeTheorem": active_range_theorem,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjoint-json", default="adjoint_eval_representer_certificate.json")
    parser.add_argument("--boundary-json", default="boundary_row_representer_certificate.json")
    parser.add_argument("--fixed-json", default="fixed_representer_theorem_scan.json")
    parser.add_argument("--closed-green-json", default="closed_trace_hardy_green_certificate.json")
    parser.add_argument(
        "--global-observability-json",
        default="global_trace_source_observability_scan.json",
    )
    parser.add_argument("--trace-to-source-json", default="trace_to_source_kernel_profile.json")
    parser.add_argument(
        "--trace-green-json",
        default="continuum_trace_to_source_green_kernel.json",
    )
    parser.add_argument(
        "--inactive-tail-json",
        default="full_theta_source_inactive_schur_tail_certificate.json",
    )
    parser.add_argument("--json-out", default="continuum_green_representer_theorem.json")
    args = parser.parse_args()

    adjoint = load_optional(args.adjoint_json)
    boundary = load_optional(args.boundary_json)
    fixed = load_optional(args.fixed_json)
    closed_green = load_optional(args.closed_green_json)
    global_obs = load_optional(args.global_observability_json)
    trace_to_source = load_optional(args.trace_to_source_json)
    trace_green = load_optional(args.trace_green_json)
    inactive_tail = load_optional(args.inactive_tail_json)

    adj = adjoint_summary(adjoint)
    bdry = boundary_summary(boundary)
    fixed_data = fixed_summary(fixed)
    closed_green_data = closed_green_summary(closed_green)
    global_obs_data = global_observability_summary(global_obs)
    trace_to_source_data = trace_to_source_summary(trace_to_source)
    trace_green_data = trace_green_summary(trace_green)
    inactive_tail_data = inactive_tail_summary(inactive_tail)

    direct_row_definition = status(
        "direct source-row definition",
        True,
        (
            "The source rows are the two explicit Lagrange/adjoint rows "
            "ell_{u,1}=B_P[h_u,.](s0) and ell_{u,2}=(P^*h_u)(s0)ev_{s0}."
        ),
    )
    finite_direct_representers = status(
        "finite direct source-row Green representers",
        bool(
            closed_green_data
            and float(closed_green_data.get("maxRangeRelativeDefect", 1.0)) < 1e-40
        ),
        (
            "Finite Galerkin sections solve the A-Riesz range equations for "
            "the actual source rows to roundoff."
        ),
    )
    active_range = status(
        "active interval trace range inclusion",
        False,
        (
            "Finite active blocks satisfy E_active=C_sample R_global to "
            "roundoff, and the sampled density refinement is stable. The "
            "continuum trace-to-source Green kernel has not yet been proved."
        ),
    )
    inactive_tail_status = status(
        "source-inactive Hardy/Schur tail",
        bool(
            inactive_tail_data
            and inactive_tail_data.get("globalSchurTailTheoremClosed") is True
        ),
        (
            "The finite full-Phi inactive source bound is below the Schur "
            "budget, but the continuum high-block Hardy/Schur tail passage is "
            "still open."
        ),
    )
    direct_green = status(
        "continuum direct Green representers",
        False,
        (
            "By Riesz, g_{u,k} exists once ell_{u,k} is A-bounded on the "
            "completed high closed-trace space. The remaining proof is active "
            "trace range inclusion plus the source-inactive tail estimate."
        ),
    )
    uniform_gram = status(
        "uniform 2x2 source Gram bound",
        False,
        (
            "Equivalent to a uniform source-row Hardy/Green bound "
            "sup_u ||E_u||_{A->C^2}<infinity. It follows after the direct "
            "representer theorem is closed."
        ),
    )
    fixed_jet_route = status(
        "fixed-jet sufficient route",
        False,
        (
            "The decomposition through k_j^hi is algebraically sufficient, "
            "but finite fixed-jet norms grow sharply; the proof should keep "
            "the special source rows together."
        ),
    )

    theorem_closed = direct_green["closed"] and uniform_gram["closed"]
    data = {
        "theoremName": "continuum Green representer construction",
        "adjointJson": args.adjoint_json if adjoint else None,
        "boundaryJson": args.boundary_json if boundary else None,
        "fixedJson": args.fixed_json if fixed else None,
        "closedGreenJson": args.closed_green_json if closed_green else None,
        "globalObservabilityJson": args.global_observability_json if global_obs else None,
        "traceToSourceJson": args.trace_to_source_json if trace_to_source else None,
        "traceGreenJson": args.trace_green_json if trace_green else None,
        "inactiveTailJson": args.inactive_tail_json if inactive_tail else None,
        "directSourceRowDefinitionStatus": direct_row_definition,
        "finiteDirectRepresenterStatus": finite_direct_representers,
        "activeTraceRangeInclusionStatus": active_range,
        "sourceInactiveTailStatus": inactive_tail_status,
        "directGreenRepresenterStatus": direct_green,
        "uniformSourceGramBoundStatus": uniform_gram,
        "fixedJetSufficientRouteStatus": fixed_jet_route,
        "continuumFixedJetRepresenterStatus": fixed_jet_route,
        "continuumGreenRepresentersClosed": theorem_closed,
        "adjointEvalCertificate": adj,
        "boundaryRowCertificate": bdry,
        "fixedRepresenterScan": fixed_data,
        "closedTraceHardyGreenCertificate": closed_green_data,
        "globalTraceSourceObservability": global_obs_data,
        "traceToSourceKernelProfile": trace_to_source_data,
        "continuumTraceToSourceGreenKernel": trace_green_data,
        "sourceInactiveTailCertificate": inactive_tail_data,
        "formalConstruction": [
            (
                "Let H_hi be the A-Hilbert completion of the high block after "
                "global trace closure and finite low/mid removal."
            ),
            (
                "Define ell_{u,1}(f)=B_P[h_u,f](s0) and "
                "ell_{u,2}(f)=(P^*h_u)(s0)f(s0)."
            ),
            (
                "If ell_{u,k} is A-bounded on H_hi, set "
                "g_{u,k}=J_A^{-1}ell_{u,k}, where J_A:H_hi->H_hi^* is the "
                "Riesz map induced by <.,.>_A."
            ),
            (
                "The active part is to be represented by an interval kernel "
                "K_u with P_delta E_u=integral K_u(a)Lambda_a da; it vanishes "
                "on ker R_global."
            ),
            (
                "The inactive part is bounded by the high-frequency "
                "Hardy/Schur tail. Hence the restricted source rows have "
                "Riesz representers in H_hi."
            ),
            (
                "For G_ij(u)=<g_{u,i},g_{u,j}>_A, "
                "lambda_max(G(u))=||E_u||_{A->C^2}^2; a uniform Hardy/Green "
                "bound gives sup_{0.08<=u<=0.52} lambda_max(G(u))<infinity."
            ),
        ],
        "finiteEvidence": {
            "maxFiniteEConstant": closed_green_data.get("maxFiniteEConstant"),
            "maxFiniteDEConstant": closed_green_data.get("maxFiniteDEConstant"),
            "maxDirectRangeRelativeDefect": closed_green_data.get("maxRangeRelativeDefect"),
            "maxGlobalHighFullFrac": global_obs_data.get("maxHighFullFrac"),
            "maxGlobalDerivativeHighFullFrac": global_obs_data.get(
                "maxDerivativeHighFullFrac"
            ),
            "traceToSourceActiveDim": trace_to_source_data.get("activeDim"),
            "traceToSourceRangeResidualRelative": trace_to_source_data.get(
                "rangeResidualRelative"
            ),
            "traceToSourceCoefficientOperator": trace_to_source_data.get(
                "coefficientOperator"
            ),
            "traceGreenSampledMaxDensityL2": trace_green_data.get("sampledMaxDensityL2"),
            "traceGreenRefinementMaxWeightedDensityL2": trace_green_data.get(
                "refinementMaxWeightedDensityL2"
            ),
            "traceGreenDensityL2Spread": trace_green_data.get(
                "refinementDensityL2RelativeSpread"
            ),
            "traceGreenSampledIbpDefect": (
                (trace_green_data.get("adjointBoundaryDiagnostic") or {}).get(
                    "maxIbpRelativeDefectOnActive"
                )
            ),
            "traceGreenRegularizedIbpDefect": (
                (trace_green_data.get("regularizedBvpDiagnostic") or {}).get(
                    "regularizedMaxIbpRelativeDefectOnActive"
                )
            ),
            "traceGreenRegularizedRangeResidual": (
                (trace_green_data.get("regularizedBvpDiagnostic") or {}).get(
                    "regularizedRangeResidualRelative"
                )
            ),
            "traceGreenJumpSolveResidual": (
                (trace_green_data.get("jumpConditionCertificate") or {}).get(
                    "maxJumpSolveResidual"
                )
            ),
            "traceGreenJumpInverseFrobenius": (
                (trace_green_data.get("jumpConditionCertificate") or {}).get(
                    "inverseJumpMatrixFrobenius"
                )
            ),
            "traceGreenInteriorSourceSolved": (
                (trace_green_data.get("jumpConditionCertificate") or {}).get(
                    "interiorSourceSolved"
                )
            ),
            "traceGreenEndpointResidual": (
                (trace_green_data.get("endpointSelectionCertificate") or {}).get(
                    "maxEndpointResidualNorm"
                )
            ),
            "traceGreenEndpointMaxZNorm": (
                (trace_green_data.get("endpointSelectionCertificate") or {}).get(
                    "maxSelectedZNorm"
                )
            ),
            "traceGreenEndpointActualRhsInRange": (
                (trace_green_data.get("endpointSelectionCertificate") or {})
                .get("actualEndpointRhsInRange", {})
                .get("closed")
            ),
            "traceGreenExactEndpointReduction": (
                (trace_green_data.get("endpointRangeTheorem") or {}).get(
                    "exactEndpointReductionClosed"
                )
            ),
            "traceGreenEndpointCompactBound": (
                (trace_green_data.get("endpointRangeTheorem") or {})
                .get("uniformTraceDualBoundFromCompactnessStatus", {})
                .get("closed")
            ),
            "traceGreenActiveAnnihilatorCriterion": (
                (trace_green_data.get("activeTraceRangeTheorem") or {})
                .get("hilbertAnnihilatorRangeCriterionStatus", {})
                .get("closed")
            ),
            "traceGreenActiveUniqueContinuation": (
                (trace_green_data.get("activeTraceRangeTheorem") or {})
                .get("closedTraceActiveUniqueContinuationStatus", {})
                .get("closed")
            ),
            "inactiveFracUpper": inactive_tail_data.get("worstContinuumInactiveFracUpper"),
            "activeTraceKernelSourceFrac": inactive_tail_data.get(
                "worstActiveTraceKernelSourceFrac"
            ),
            "evalRepresenterNorm2Basis20": adj.get("evalRepresenterNorm2"),
            "boundaryCoverNorm2Basis20": bdry.get("coverBoundaryNorm2"),
            "evalCoverNorm2Basis20": adj.get("coverEvalNorm2"),
            "fixedScanMaxEvalRepresenterNorm2": fixed_data.get("maxEvalRepresenterNorm2"),
            "fixedScanMaxBoundaryCoverNorm2": fixed_data.get("maxBoundaryCoverNorm2"),
            "maxFactorRelativeDefect": max(
                float(adj.get("maxFactorRelativeDefect", 0.0) or 0.0),
                float(bdry.get("maxFactorRelativeDefect", 0.0) or 0.0),
                float(fixed_data.get("maxJetRangeRelativeDefect", 0.0) or 0.0),
            ),
        },
        "remainingAnalyticGap": (
            "Prove the continuum active trace-to-source range inclusion and "
            "the continuum source-inactive Hardy/Schur tail bound for the "
            "actual source rows. Do not require separate boundedness of all "
            "endpoint jet functionals."
        ),
        "interpretation": (
            "The continuum Green representers should be constructed directly "
            "as A-Riesz representers of the two Lagrange source rows. Finite "
            "direct representers and active trace-to-source factorizations are "
            "exact to roundoff; the continuum range/tail theorem remains open."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum Green representer construction")
    print(f"  direct row definition closed: {direct_row_definition['closed']}")
    print(f"  finite direct representers: {finite_direct_representers['closed']}")
    print(f"  active trace range inclusion: {active_range['closed']}")
    print(f"  inactive tail theorem: {inactive_tail_status['closed']}")
    if closed_green_data:
        print(
            "  finite max ||E_u||_A^2: "
            f"{float(closed_green_data['maxFiniteEConstant']):.12e}"
        )
    if global_obs_data:
        print(
            "  max global high/full frac: "
            f"{float(global_obs_data['maxHighFullFrac']):.12e}"
        )
    if trace_to_source_data:
        print(
            "  trace-to-source rel residual: "
            f"{float(trace_to_source_data['rangeResidualRelative']):.12e}"
        )
    if trace_green_data:
        print(
            "  trace Green density L2 spread: "
            f"{float(trace_green_data['refinementDensityL2RelativeSpread']):.12e}"
        )
    if adj:
        print(f"  eval ||k0||_A^2 basis20: {float(adj['evalRepresenterNorm2']):.12e}")
    if fixed_data:
        print(
            "  fixed scan max ||k0||_A^2: "
            f"{float(fixed_data['maxEvalRepresenterNorm2']):.12e}"
        )
    print(f"  continuum Green representers closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
