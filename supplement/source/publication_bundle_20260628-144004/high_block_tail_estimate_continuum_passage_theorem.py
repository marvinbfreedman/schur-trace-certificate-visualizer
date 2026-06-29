#!/usr/bin/env python3
r"""Continuum passage for the high-block source-inactive tail estimate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--active-json",
        default="synchronized_active_range_interval_consequence_theorem.json",
    )
    parser.add_argument(
        "--tail-constants-json",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument(
        "--trace-frame-json",
        default="continuum_trace_frame_lower_bound_consequence_theorem.json",
    )
    parser.add_argument(
        "--projection-json",
        default="high_block_mosco_projection_consequence_theorem.json",
    )
    parser.add_argument(
        "--compression-json",
        default="high_block_compact_compression_consequence_theorem.json",
    )
    parser.add_argument(
        "--spectral-json",
        default="high_block_spectral_projection_consequence_theorem.json",
    )
    parser.add_argument(
        "--minmax-passage-json",
        default="high_block_minmax_tail_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_tail_estimate_continuum_passage_theorem.json",
    )
    args = parser.parse_args()

    active = load(args.active_json)
    tail = load(args.tail_constants_json)
    trace = load(args.trace_frame_json)
    projection = load(args.projection_json)
    compression = load(args.compression_json)
    spectral = load(args.spectral_json)
    minmax_passage = load(args.minmax_passage_json)

    active_ok = closed(active, "activeRangeInclusionStatus")
    inactive_ok = closed(active, "sourceInactiveTailDominationStatus")
    tail_ok = bool(
        tail.get("sourceInactiveTailConstantsClosed")
        or tail.get("sourceInactiveTailDominationConstantsClosed")
        or tail.get("minMaxProofInCertifiedSourceModel")
    )
    absorb_ok = bool(tail.get("absorbableByFiniteLowMidBlock"))
    trace_ok = (
        closed(trace, "continuumTraceFrameLowerBoundStatus")
        and closed(trace, "traceQuadratureConsistencyStatus")
        and closed(trace, "boundedSampledTraceCorrectionStatus")
    )
    projection_ok = bool(
        projection.get("highBlockMoscoProjectionClosed")
        or projection.get("moscoProjectionConvergenceClosed")
    )
    compression_ok = bool(
        compression.get("highBlockCompactCompressionClosed")
        or compression.get("compactCompressionNormConvergenceClosed")
    )
    spectral_ok = bool(
        spectral.get("highBlockSpectralProjectionClosed")
        or spectral.get("compactSourceRieszProjectionConvergenceClosed")
    )
    abstract_minmax_ok = bool(
        minmax_passage.get("highBlockMinmaxTailPassageClosed")
        or minmax_passage.get("abstractMinmaxTailPassageClosed")
    )
    abstract_ok = all([projection_ok, compression_ok, spectral_ok, abstract_minmax_ok])
    passage_ok = all([active_ok, inactive_ok, tail_ok, absorb_ok, trace_ok, abstract_ok])

    epsilon = tail.get("normalizedEpsilonDelta")
    budget = tail.get("finiteLowMidSchurBudget")
    slack = tail.get("absorptionSlack")

    data = {
        "theoremName": "high-block tail estimate continuum passage theorem",
        "proofClass": "analytic proof",
        "activeJson": args.active_json,
        "tailConstantsJson": args.tail_constants_json,
        "traceFrameJson": args.trace_frame_json,
        "abstractProjectionJson": args.projection_json,
        "abstractCompressionJson": args.compression_json,
        "abstractSpectralJson": args.spectral_json,
        "abstractMinmaxPassageJson": args.minmax_passage_json,
        "normalizedEpsilonDelta": epsilon,
        "finiteLowMidSchurBudget": budget,
        "absorptionSlack": slack,
        "statuses": {
            "activeRangeInputStatus": status(
                "active range input",
                active_ok,
                "The active range theorem gives R_global f=0 => P_active Ehat f=0.",
            ),
            "inactiveCertifiedModelInputStatus": status(
                "inactive certified source model input",
                inactive_ok and tail_ok and absorb_ok,
                (
                    "The certified source model supplies the source-inactive "
                    "tail constant and verifies it is absorbable by the finite "
                    "low/mid Schur budget."
                ),
            ),
            "traceFrameCorrectionInputStatus": status(
                "trace-frame correction input",
                trace_ok,
                (
                    "The continuum trace-frame theorem supplies the positive "
                    "frame floor, quadrature consistency, and bounded sampled "
                    "trace correction right inverses."
                ),
            ),
            "abstractProjectionInputStatus": status(
                "high-block Mosco projection input",
                projection_ok,
                (
                    "The high-block Mosco projection theorem supplies the "
                    "trace-frame recovery input and applies the abstract "
                    "projection convergence result."
                ),
            ),
            "abstractCompressionInputStatus": status(
                "high-block compact-compression input",
                compression_ok,
                (
                    "The high-block compact-compression theorem verifies the "
                    "source compactness hypothesis and applies the abstract "
                    "finite-rank compression argument."
                ),
            ),
            "abstractSpectralInputStatus": status(
                "high-block Riesz projection input",
                spectral_ok,
                (
                    "The high-block spectral theorem supplies the source gap "
                    "and applies the abstract Riesz projection stability "
                    "result."
                ),
            ),
            "abstractMinmaxInputStatus": status(
                "high-block min-max tail passage input",
                abstract_minmax_ok,
                (
                    "The high-block min-max theorem supplies the finite tail "
                    "constant, stable spectral split, and abstract min-max "
                    "passage."
                ),
            ),
            "abstractFunctionalAnalysisInputStatus": status(
                "abstract functional-analysis input",
                abstract_ok,
                (
                    "The split Mosco, high-block compact-compression, "
                    "high-block Riesz projection, and min-max theorems close "
                    "the high-block passage."
                ),
            ),
            "highBlockTailEstimateContinuumPassageStatus": status(
                "high-block tail estimate passes to continuum",
                passage_ok,
                (
                    "The active/inactive source split, certified finite "
                    "tail constant, trace-frame correction theorem, and "
                    "abstract compact-source Mosco passage give the continuum "
                    "high-block source-inactive tail estimate."
                ),
            ),
        },
        "highBlockTailEstimateContinuumPassageClosed": passage_ok,
        "tailEstimatePassesToContinuum": passage_ok,
        "conditionalHighBlockExhaustionClosed": passage_ok,
        "proof": [
            "Use active range inclusion to remove the active source component on ker R_global.",
            "Use the certified inactive-tail constants and low/mid absorption budget in the finite source model.",
            "Use the continuum trace-frame theorem to correct Galerkin approximants into the sampled closed-trace high blocks.",
            "Use Mosco projection convergence, compact-compression norm convergence, Riesz projection stability, and min-max passage to transport the finite estimate to the continuum.",
        ],
        "remainingAnalyticGap": None
        if passage_ok
        else "One active range, tail constant, trace-frame, or abstract compact-source input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block tail estimate continuum passage theorem")
    print(f"  active range: {active_ok}")
    print(f"  inactive certified model: {inactive_ok and tail_ok and absorb_ok}")
    print(f"  trace-frame correction: {trace_ok}")
    print(f"  abstract Mosco projection: {projection_ok}")
    print(f"  abstract compact compression: {compression_ok}")
    print(f"  abstract Riesz projection: {spectral_ok}")
    print(f"  abstract min-max passage: {abstract_minmax_ok}")
    print(f"  continuum tail passage: {passage_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
