#!/usr/bin/env python3
r"""Assemble the direct and closed-cone KLM-to-de Branges bridge attempts."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hardy-json", default="klm_debranges_canonical_hardy_image.json")
    parser.add_argument("--pullback-json", default="klm_debranges_pullback_probe.json")
    parser.add_argument("--green-json", default="continuum_green_lift_closure_theorem.json")
    parser.add_argument("--branch-json", default="klm_debranges_branch_transport_theorem.json")
    parser.add_argument("--trace-map-json", default="klm_debranges_trace_map_constructor.json")
    parser.add_argument("--lambda-json", default="klm_debranges_lambda_trace_candidate.json")
    parser.add_argument("--feature-inverse-json", default="klm_debranges_feature_inverse_candidate.json")
    parser.add_argument("--continuous-normal-json", default="klm_debranges_continuous_normal_equation.json")
    parser.add_argument("--coupled-normal-json", default="klm_debranges_coupled_branch_normal_equation.json")
    parser.add_argument("--transmutation-json", default="klm_debranges_transmutation_kernel_probe.json")
    parser.add_argument("--nonlocal-transmutation-json", default="klm_debranges_nonlocal_transmutation_probe.json")
    parser.add_argument("--mode-match-json", default="xi_mellin_volterra_mode_match.json")
    parser.add_argument("--mode-mixing-json", default="xi_mellin_volterra_mode_mixing.json")
    parser.add_argument("--mellin-boundary-json", default="xi_mellin_convolution_boundary_identity.json")
    parser.add_argument("--boundary-trace-json", default="xi_boundary_prefix_trace_resolution.json")
    parser.add_argument("--boundary-concomitant-json", default="xi_mellin_boundary_concomitant.json")
    parser.add_argument("--augmented-repair-json", default="xi_augmented_trace_repair_schur.json")
    parser.add_argument("--augmented-continuum-json", default="xi_augmented_trace_continuum_lift.json")
    parser.add_argument(
        "--augmented-closed-cone-json",
        default="klm_debranges_augmented_closed_cone_theorem.json",
    )
    parser.add_argument("--json-out", default="klm_debranges_bridge_attempt.json")
    args = parser.parse_args()

    hardy = load(args.hardy_json)
    pullback = load(args.pullback_json)
    green = load(args.green_json)
    branch = load(args.branch_json) if Path(args.branch_json).exists() else {}
    trace_map = load(args.trace_map_json) if Path(args.trace_map_json).exists() else {}
    lambda_candidate = load(args.lambda_json) if Path(args.lambda_json).exists() else {}
    feature_inverse = load(args.feature_inverse_json) if Path(args.feature_inverse_json).exists() else {}
    continuous_normal = load(args.continuous_normal_json) if Path(args.continuous_normal_json).exists() else {}
    coupled_normal = load(args.coupled_normal_json) if Path(args.coupled_normal_json).exists() else {}
    transmutation = load(args.transmutation_json) if Path(args.transmutation_json).exists() else {}
    nonlocal_transmutation = (
        load(args.nonlocal_transmutation_json) if Path(args.nonlocal_transmutation_json).exists() else {}
    )
    mode_match = load(args.mode_match_json) if Path(args.mode_match_json).exists() else {}
    mode_mixing = load(args.mode_mixing_json) if Path(args.mode_mixing_json).exists() else {}
    mellin_boundary = load(args.mellin_boundary_json) if Path(args.mellin_boundary_json).exists() else {}
    boundary_trace = load(args.boundary_trace_json) if Path(args.boundary_trace_json).exists() else {}
    boundary_concomitant = (
        load(args.boundary_concomitant_json) if Path(args.boundary_concomitant_json).exists() else {}
    )
    augmented_repair = load(args.augmented_repair_json) if Path(args.augmented_repair_json).exists() else {}
    augmented_continuum = (
        load(args.augmented_continuum_json) if Path(args.augmented_continuum_json).exists() else {}
    )
    augmented_closed_cone = (
        load(args.augmented_closed_cone_json)
        if Path(args.augmented_closed_cone_json).exists()
        else {}
    )

    hardy_closed = bool(hardy.get("canonicalHardyImageClosed"))
    direct_closed = bool(pullback.get("foundExactFinitePullback"))
    positive_residual = bool(pullback.get("foundPositiveResidualDomination"))
    green_closed = bool(
        green.get("greenLiftContractionClosed")
        or green.get("continuumGreenLiftClosureClosed")
        or green.get("statuses", {}).get("greenLiftContractionStatus", {}).get("closed")
    )
    contraction = hardy.get("hardyBranchContraction", {})

    statuses = {
        "canonicalHardyImageStatus": status(
            "canonical de Branges Hardy image",
            hardy_closed,
            (
                "The shifted-Xi kernel is exactly the signed Hardy branch Gram "
                "K_E=<h^+,h^+>-<h^-,h^->.  The numerical implementation matches "
                f"the direct kernel with relative error {hardy.get('hardyExactRelativeError')}."
            ),
            blocker=None if hardy_closed else "Fix the shifted-Xi/Hardy normalization.",
        ),
        "directKlmPullbackStatus": status(
            "direct finite coherent KLM pullback",
            direct_closed,
            (
                "The tested Gaussian/coherent phase-space ansatzes did not produce "
                "an exact pullback.  The best residual was indefinite, so the direct "
                "naive coherent embedding is rejected."
            ),
            blocker="Do not tune Gaussian packets further; derive the canonical Weyl image.",
        ),
        "positiveResidualDominationStatus": status(
            "finite positive residual domination",
            positive_residual,
            (
                "The finite probe did not find K_E - T^*KLM T >= 0 for the tested "
                "coherent ansatzes."
            ),
            blocker="Any domination proof needs the canonical branch transport, not the tested packets.",
        ),
        "volterraGreenContractionInputStatus": status(
            "completed Volterra/KLM branch contraction input",
            green_closed,
            (
                "The earlier completed Green-lift theorem gives the Volterra branch "
                "contraction G_- = C K E G_+ on the completed trace-fiber domain."
            ),
            blocker=None if green_closed else "Close continuum_green_lift_closure_theorem.py.",
        ),
        "closedConeBridgeStatus": status(
            "closed positive-cone KLM-to-de Branges bridge",
            bool(augmented_closed_cone.get("bridgeClosed")),
            (
                "The augmented closed-cone theorem now supplies the missing "
                "evaluation pullback in closed-cone form using "
                "R_aug=(Lambda,Mu).  The exact Mellin split gives the finite "
                "theta pullbacks, the augmented continuum repair puts those "
                "pullbacks in the KLM positive cone, and compact theta-tail "
                "convergence identifies the limit with the shifted-Xi de "
                "Branges kernel.  Historically, this route was reduced through "
                "the exact-U criterion and Hardy strong truncation limit; before "
                "the augmented trace was found, what remained was "
                "the concrete evaluation trace map z -> x_z, or truncated maps "
                "x_{z,R}, whose Volterra/KLM feature Grams converge to the Hardy "
                "branch Grams.  The first finite constructor shows that abstract "
                "plus-Gram lifting and naive profile fitting do not supply this map; "
                "the Lambda trace candidate shows that bare exponential endpoint "
                "jets do not supply it either.  The coefficient-level finite "
                "inverse of G_+ was also tested and does not close the map in "
                "the sampled Weyl/Volterra feature model.  The continuous "
                "Galerkin normal equation G_+^*G_+f_z=G_+^*h_z^+ was then "
                "assembled from the exact Weyl/Volterra kernel and solved to "
                "roundoff; its Lambda traces still do not recover the joint "
                "Hardy branch Grams.  The coupled signed-branch equation was "
                "also solved for direct and swapped branch normalizations, "
                "again without recovering the Hardy joint Gram.  The natural "
                "coarea transmutation ansatz r=s+u was then tested with scalar "
                "weights theta(s,u,sigma); it improves neither the lifted nor "
                "the trace-side bridge enough to close the map.  A finite "
                "nonlocal resolvent/heat dictionary centered at the natural "
                "Volterra/Mellin locations gives only a modest improvement and "
                "still does not close the bridge.  The exact finite-core "
                "incomplete-gamma Xi atoms were then matched term-by-term with "
                "the exact Volterra ratio atoms; the Xi atom normalization is "
                "correct, but the diagonal atom dictionary also fails to close "
                "the joint Gram.  A finite 6x6 least-squares mode-mixing matrix "
                "between the incomplete-gamma Xi atoms and the Volterra ratio "
                "atoms improves the lifted residual only modestly and is badly "
                "conditioned, so it is a numerical projection rather than a "
                "structural intertwiner.  The exact Mellin split explains the "
                "failure: the true identity is a diagonal moving Volterra tail "
                "plus an s-dependent incomplete-gamma boundary prefix, not a "
                "constant atom-mixing matrix.  A finite trace-resolution test "
                "shows that this boundary prefix is not killed by the current "
                "sampled Lambda_a endpoint trace nullspace.  The missing "
                "Mellin-boundary concomitant is now identified as the Volterra "
                "primitive trace Mu_z(f)=int B(s,z)f(s)ds; augmenting Lambda_a "
                "by Mu_z kills the boundary prefix on the finite augmented "
                "trace nullspace.  The finite augmented Schur repair is positive "
                "and annihilates Mu on ker(Lambda,Mu).  The continuum lift of "
                "this augmented repair is now closed on the completed augmented "
                "Volterra/Mellin trace-fiber domain, and the final augmented "
                "closed-cone theorem closes the pullback limit."
            ),
            blocker=None
            if bool(augmented_closed_cone.get("bridgeClosed"))
            else (
                augmented_continuum.get("nextProofTarget")
                or augmented_repair.get("nextProofTarget")
                or boundary_concomitant.get("nextProofTarget")
                or boundary_trace.get("nextProofTarget")
                or mellin_boundary.get("nextProofTarget")
                or mode_mixing.get("nextProofTarget")
                or mode_match.get("nextProofTarget")
                or nonlocal_transmutation.get("nextProofTarget")
                or transmutation.get("nextProofTarget")
                or coupled_normal.get("nextProofTarget")
                or continuous_normal.get("nextProofTarget")
                or feature_inverse.get("nextProofTarget")
                or "Derive the missing Hardy-to-Volterra transform or branch normalization."
            ),
        ),
    }

    data = {
        "theoremName": "direct and closed-cone KLM-to-de Branges bridge attempt",
        "hardyJson": args.hardy_json,
        "pullbackJson": args.pullback_json,
        "greenJson": args.green_json,
        "branchJson": args.branch_json,
        "traceMapJson": args.trace_map_json,
        "lambdaTraceJson": args.lambda_json,
        "featureInverseJson": args.feature_inverse_json,
        "continuousNormalJson": args.continuous_normal_json,
        "coupledNormalJson": args.coupled_normal_json,
        "transmutationJson": args.transmutation_json,
        "nonlocalTransmutationJson": args.nonlocal_transmutation_json,
        "modeMatchJson": args.mode_match_json,
        "modeMixingJson": args.mode_mixing_json,
        "mellinBoundaryJson": args.mellin_boundary_json,
        "boundaryTraceJson": args.boundary_trace_json,
        "boundaryConcomitantJson": args.boundary_concomitant_json,
        "augmentedRepairJson": args.augmented_repair_json,
        "augmentedContinuumJson": args.augmented_continuum_json,
        "augmentedClosedConeJson": args.augmented_closed_cone_json,
        "statuses": statuses,
        "hardyBranchContractionSample": contraction,
        "branchTransportTheorem": {
            "available": bool(branch),
            "exactUCriterionClosed": bool(
                branch.get("statuses", {}).get("exactUCriterionStatus", {}).get("closed")
            ),
            "exactUConstructed": bool(
                branch.get("statuses", {}).get("exactUConstructedStatus", {}).get("closed")
            ),
            "hardyStrongTruncationClosed": bool(
                branch.get("statuses", {}).get("hardyStrongTruncationStatus", {}).get("closed")
            ),
            "klmClosedConeLimitClosed": bool(
                branch.get("statuses", {}).get("klmClosedConeLimitStatus", {}).get("closed")
            ),
            "hardyStrongLimitProof": branch.get("hardyStrongLimitProof"),
        },
        "traceMapConstructor": {
            "available": bool(trace_map),
            "exactTraceMapConstructed": bool(trace_map.get("exactTraceMapConstructed")),
            "truncatedTraceMapConstructed": bool(trace_map.get("truncatedTraceMapConstructed")),
            "fullHardyTarget": trace_map.get("fullHardyTarget"),
            "truncatedHardyTarget": trace_map.get("truncatedHardyTarget"),
            "profileFitCandidates": trace_map.get("profileFitCandidates"),
            "diagnosis": trace_map.get("diagnosis"),
        },
        "lambdaTraceCandidate": {
            "available": bool(lambda_candidate),
            "exactTraceMapConstructed": bool(lambda_candidate.get("exactTraceMapConstructed")),
            "bestByScaledFullJointResidual": lambda_candidate.get("bestByScaledFullJointResidual"),
            "diagnosis": lambda_candidate.get("diagnosis"),
        },
        "featureInverseCandidate": {
            "available": bool(feature_inverse),
            "directFeatureMapConstructed": bool(feature_inverse.get("directFeatureMapConstructed")),
            "traceMapConstructed": bool(feature_inverse.get("traceMapConstructed")),
            "bestDirectFeatureGram": feature_inverse.get("bestDirectFeatureGram"),
            "bestGreenTraceGram": feature_inverse.get("bestGreenTraceGram"),
            "diagnosis": feature_inverse.get("diagnosis"),
        },
        "continuousNormalEquation": {
            "available": bool(continuous_normal),
            "continuousNormalEquationSolved": bool(continuous_normal.get("continuousNormalEquationSolved")),
            "directContinuousMapConstructed": bool(continuous_normal.get("directContinuousMapConstructed")),
            "traceMapConstructed": bool(continuous_normal.get("traceMapConstructed")),
            "bestDirectContinuousGram": continuous_normal.get("bestDirectContinuousGram"),
            "bestGreenTraceGram": continuous_normal.get("bestGreenTraceGram"),
            "diagnosis": continuous_normal.get("diagnosis"),
        },
        "coupledBranchNormalEquation": {
            "available": bool(coupled_normal),
            "coupledNormalEquationSolved": bool(coupled_normal.get("coupledNormalEquationSolved")),
            "directContinuousMapConstructed": bool(coupled_normal.get("directContinuousMapConstructed")),
            "traceMapConstructed": bool(coupled_normal.get("traceMapConstructed")),
            "bestDirectContinuousGram": coupled_normal.get("bestDirectContinuousGram"),
            "bestGreenTraceGram": coupled_normal.get("bestGreenTraceGram"),
            "diagnosis": coupled_normal.get("diagnosis"),
        },
        "transmutationKernelProbe": {
            "available": bool(transmutation),
            "simpleCoareaTransmutationConstructed": bool(
                transmutation.get("simpleCoareaTransmutationConstructed")
            ),
            "bestIntegratedVolterraGram": transmutation.get("bestIntegratedVolterraGram"),
            "bestLiftedTriangleGram": transmutation.get("bestLiftedTriangleGram"),
            "bestGreenTraceGram": transmutation.get("bestGreenTraceGram"),
            "diagnosis": transmutation.get("diagnosis"),
        },
        "nonlocalTransmutationKernelProbe": {
            "available": bool(nonlocal_transmutation),
            "nonlocalDictionaryConstructedBridge": bool(
                nonlocal_transmutation.get("nonlocalDictionaryConstructedBridge")
            ),
            "bestLiftedTriangleGram": nonlocal_transmutation.get("bestLiftedTriangleGram"),
            "bestIntegratedVolterraGram": nonlocal_transmutation.get("bestIntegratedVolterraGram"),
            "bestGreenTraceGram": nonlocal_transmutation.get("bestGreenTraceGram"),
            "diagnosis": nonlocal_transmutation.get("diagnosis"),
        },
        "xiMellinVolterraModeMatch": {
            "available": bool(mode_match),
            "modeMatchedBridgeConstructed": bool(mode_match.get("modeMatchedBridgeConstructed")),
            "maxXiAtomDiagnosticRelativeError": mode_match.get("maxXiAtomDiagnosticRelativeError"),
            "bestLiftedTriangleGram": mode_match.get("bestLiftedTriangleGram"),
            "bestIntegratedVolterraGram": mode_match.get("bestIntegratedVolterraGram"),
            "bestGreenTraceGram": mode_match.get("bestGreenTraceGram"),
            "diagnosis": mode_match.get("diagnosis"),
        },
        "xiMellinVolterraModeMixing": {
            "available": bool(mode_mixing),
            "modeMixingBridgeConstructed": bool(mode_mixing.get("modeMixingBridgeConstructed")),
            "maxXiAtomDiagnosticRelativeError": mode_mixing.get("maxXiAtomDiagnosticRelativeError"),
            "bestMixingFit": mode_mixing.get("bestMixingFit"),
            "bestLiftedTriangleGram": mode_mixing.get("bestLiftedTriangleGram"),
            "bestIntegratedVolterraGram": mode_mixing.get("bestIntegratedVolterraGram"),
            "bestGreenTraceGram": mode_mixing.get("bestGreenTraceGram"),
            "diagnosis": mode_mixing.get("diagnosis"),
        },
        "xiMellinConvolutionBoundaryIdentity": {
            "available": bool(mellin_boundary),
            "maxTotalSplitRelativeError": mellin_boundary.get("maxTotalSplitRelativeError"),
            "maxTailVolterraRelativeError": mellin_boundary.get("maxTailVolterraRelativeError"),
            "maxBoundaryFraction": mellin_boundary.get("maxBoundaryFraction"),
            "minTailFraction": mellin_boundary.get("minTailFraction"),
            "maxTailFraction": mellin_boundary.get("maxTailFraction"),
            "conclusion": mellin_boundary.get("conclusion"),
        },
        "xiBoundaryPrefixTraceResolution": {
            "available": bool(boundary_trace),
            "boundaryPrefixTraceResolved": bool(boundary_trace.get("boundaryPrefixTraceResolved")),
            "boundaryPrefixKilledOnTraceNullspace": bool(
                boundary_trace.get("boundaryPrefixKilledOnTraceNullspace")
            ),
            "maxRowSpanResidualRelative": boundary_trace.get("maxRowSpanResidualRelative"),
            "maxNullEnergyOpMax": boundary_trace.get("maxNullEnergyOpMax"),
            "maxNullEnergyOpL2": boundary_trace.get("maxNullEnergyOpL2"),
            "diagnosis": boundary_trace.get("diagnosis"),
        },
        "xiMellinBoundaryConcomitant": {
            "available": bool(boundary_concomitant),
            "mellinBoundaryConcomitantDerived": bool(
                boundary_concomitant.get("mellinBoundaryConcomitantDerived")
            ),
            "augmentedTraceKillsBoundaryPrefix": bool(
                boundary_concomitant.get("augmentedTraceKillsBoundaryPrefix")
            ),
            "oldTraceBoundaryNullEnergyOpMax": boundary_concomitant.get(
                "oldTraceBoundaryNullEnergyOpMax"
            ),
            "augmentedTraceBoundaryNullEnergyOpMax": boundary_concomitant.get(
                "augmentedTraceBoundaryNullEnergyOpMax"
            ),
            "concomitantIdentity": boundary_concomitant.get("concomitantIdentity"),
            "diagnosis": boundary_concomitant.get("diagnosis"),
        },
        "xiAugmentedTraceRepairSchur": {
            "available": bool(augmented_repair),
            "augmentedRepairPositive": bool(augmented_repair.get("augmentedRepairPositive")),
            "augmentedMuAnnihilated": bool(augmented_repair.get("augmentedMuAnnihilated")),
            "oldMuActionOnLambdaNullspace": augmented_repair.get("oldMuActionOnLambdaNullspace"),
            "augmentedMuActionOnAugmentedNullspace": augmented_repair.get(
                "augmentedMuActionOnAugmentedNullspace"
            ),
            "augmentedRepair": augmented_repair.get("augmentedRepair"),
            "diagnosis": augmented_repair.get("diagnosis"),
        },
        "xiAugmentedTraceContinuumLift": {
            "available": bool(augmented_continuum),
            "continuumAugmentedRepairClosed": bool(
                augmented_continuum.get("continuumAugmentedRepairClosed")
            ),
            "statuses": augmented_continuum.get("statuses"),
            "scope": augmented_continuum.get("scope"),
        },
        "augmentedClosedConeTheorem": {
            "available": bool(augmented_closed_cone),
            "bridgeClosed": bool(augmented_closed_cone.get("bridgeClosed")),
            "transformClosed": bool(augmented_closed_cone.get("transformClosed")),
            "sampleConstants": augmented_closed_cone.get("sampleConstants"),
            "nextProofTarget": augmented_closed_cone.get("nextProofTarget"),
        },
        "bestFinitePullbackCandidate": pullback.get("bestCandidate"),
        "directKlmPullbackClosed": direct_closed,
        "closedPositiveConeBridgeClosed": bool(augmented_closed_cone.get("bridgeClosed")),
        "rhClosed": False,
        "nextProofTarget": (
            augmented_closed_cone.get("nextProofTarget")
            if bool(augmented_closed_cone.get("bridgeClosed"))
            else statuses["closedConeBridgeStatus"]["blocker"]
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("KLM-to-de Branges bridge attempt")
    print(f"  canonical Hardy image closed: {hardy_closed}")
    print(f"  direct coherent KLM pullback closed: {direct_closed}")
    print(f"  positive residual domination found: {positive_residual}")
    print(f"  Volterra/KLM contraction input closed: {green_closed}")
    if contraction.get("available"):
        print(f"  sample Hardy contraction max eigenvalue: {contraction['maxGeneralizedEigenvalue']:.12e}")
        print(f"  sample Hardy contraction margin: {contraction['sampleMargin']:.12e}")
    print(f"  closed-cone bridge closed: {bool(augmented_closed_cone.get('bridgeClosed'))}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
