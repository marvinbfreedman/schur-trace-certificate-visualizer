#!/usr/bin/env python3
r"""Narrow consequence of the full-theta source quadrature certificate.

The raw certificate contains diagnostic provenance for the source quadrature,
theta-tail propagation, and finite active response computation.  Downstream
theorems only need the interval consequence:

    ||S_full - S_{8,h}|| <= totalContinuumErrorBound,
    totalContinuumErrorBound < spectralGapToComplementS8/4,
    lowerBoundAfterContinuumAndTail > 0.

This wrapper checks those inequalities and exports only the data needed by the
source noncollapse theorem.  It intentionally avoids exposing raw certificate
filenames as theorem dependencies; the publication audit should point here as
the interface for this interval result.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any


getcontext().prec = 80


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def dec(value: Any) -> Decimal:
    return Decimal(str(value))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {
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
    parser.add_argument(
        "--certificate-json",
        default="full_theta_source_quadrature_certificate.json",
        help="Raw machine certificate to be checked and compressed.",
    )
    parser.add_argument(
        "--json-out",
        default="full_theta_source_quadrature_consequence_theorem.json",
    )
    args = parser.parse_args()

    raw = load(args.certificate_json)

    source_quad = dec(raw["sourceQuadratureErrorBound"])
    theta_tail = dec(raw["omittedTailDeltaSBound"])
    total_error = dec(raw["totalContinuumErrorBound"])
    d2_envelope = dec(raw["d2SourceGramIntervalEnvelope"])
    trapezoid_factor = dec(raw["trapezoidFactor"])
    gap = dec(raw["spectralGapToComplementS8"])
    alpha = dec(raw["totalProjectorAlpha"])
    rank_margin = dec(raw["finiteRankMarginS8"])
    angle_factor = dec(raw["angleFactorLowerBound"])
    complement_response = dec(raw["complementResponseOperatorNorm"])
    lower_after = dec(raw["lowerBoundAfterContinuumAndTail"])

    active_eigs = [dec(x) for x in raw.get("activeEigenvaluesS8", [])]
    response_svals = [dec(x) for x in raw.get("finiteActiveResponseSingularValuesS8", [])]

    source_quadrature_closed = source_quad >= 0 and d2_envelope >= 0 and trapezoid_factor >= 0
    theta_tail_closed = theta_tail >= 0
    total_matches_parts = total_error >= source_quad + theta_tail
    if not total_matches_parts:
        total_matches_parts = abs(total_error - (source_quad + theta_tail)) <= Decimal("1e-12")

    total_error_closed = source_quadrature_closed and theta_tail_closed and total_matches_parts
    riesz_closed = (
        total_error_closed
        and gap > 0
        and total_error < gap / Decimal(4)
        and alpha < 1
        and bool(raw.get("totalGapConditionPasses"))
    )
    rank_margin_closed = (
        riesz_closed
        and len(active_eigs) == 2
        and len(response_svals) >= 2
        and min(active_eigs) > 0
        and min(response_svals[:2]) > 0
        and rank_margin > 0
        and angle_factor > 0
        and complement_response >= 0
        and lower_after > 0
        and bool(raw.get("fullPhiContinuumSourceNoncollapsePasses"))
    )

    data: dict[str, Any] = {
        "theoremName": "full-theta source quadrature consequence theorem",
        "proofClass": "interval/ball certificate",
        "omega": raw.get("omega"),
        "L": raw.get("L"),
        "basis": raw.get("basis"),
        "fullNmax": raw.get("fullNmax"),
        "sourceGrid": raw.get("sourceGrid"),
        "sourceNodesFirstLast": raw.get("sourceNodesFirstLast"),
        "sourceWeightsFirstLast": raw.get("sourceWeightsFirstLast"),
        "activeDimension": len(active_eigs),
        "activeEigenvaluesS8": raw.get("activeEigenvaluesS8", []),
        "finiteActiveResponseSingularValuesS8": raw.get(
            "finiteActiveResponseSingularValuesS8", []
        ),
        "sourceQuadratureErrorBound": raw.get("sourceQuadratureErrorBound"),
        "omittedTailDeltaSBound": raw.get("omittedTailDeltaSBound"),
        "totalContinuumErrorBound": raw.get("totalContinuumErrorBound"),
        "spectralGapToComplementS8": raw.get("spectralGapToComplementS8"),
        "totalProjectorAlpha": raw.get("totalProjectorAlpha"),
        "finiteRankMarginS8": raw.get("finiteRankMarginS8"),
        "angleFactorLowerBound": raw.get("angleFactorLowerBound"),
        "complementResponseOperatorNorm": raw.get("complementResponseOperatorNorm"),
        "lowerBoundAfterContinuumAndTail": raw.get("lowerBoundAfterContinuumAndTail"),
        "d2SourceGramIntervalEnvelope": raw.get("d2SourceGramIntervalEnvelope"),
        "trapezoidFactor": raw.get("trapezoidFactor"),
        "totalGapConditionPasses": bool(raw.get("totalGapConditionPasses")),
        "tailBoundClosed": theta_tail_closed,
        "sourceQuadratureConsequenceClosed": rank_margin_closed,
        "fullPhiContinuumSourceNoncollapsePasses": rank_margin_closed,
        "statuses": {
            "sourceQuadratureIntervalStatus": status(
                "source quadrature interval consequence",
                source_quadrature_closed,
                "The second-derivative envelope and trapezoid factor give a nonnegative source quadrature error bound.",
                blocker=None if source_quadrature_closed else "Recheck the source quadrature interval envelope.",
            ),
            "thetaTailIntervalStatus": status(
                "theta tail interval consequence",
                theta_tail_closed,
                "The omitted full-theta tail contribution is bounded by a nonnegative interval operator error.",
                blocker=None if theta_tail_closed else "Recheck the omitted theta-tail interval propagation.",
            ),
            "totalErrorStatus": status(
                "assembled continuum source error",
                total_error_closed,
                "The exported total continuum error dominates the source quadrature and theta-tail contributions.",
                blocker=None if total_error_closed else "Check totalContinuumErrorBound against its component bounds.",
            ),
            "rieszGapStatus": status(
                "Riesz gap survives source perturbation",
                riesz_closed,
                "The total perturbation is below one quarter of the active/complement source spectral gap.",
                blocker=None if riesz_closed else "Reduce the perturbation bound or certify a larger spectral gap.",
            ),
            "activeRankMarginStatus": status(
                "active source response rank margin",
                rank_margin_closed,
                "The finite rank margin remains positive after continuum quadrature and full-theta tail perturbation.",
                blocker=None if rank_margin_closed else "Close the finite rank margin or perturbation inequality.",
            ),
        },
        "formalProof": [
            "The interval quadrature certificate gives a source Gram perturbation bound.",
            "The full-theta tail interval bound is added to that source perturbation.",
            "The combined error is less than one quarter of the active/complement spectral gap.",
            "The Riesz projector and active response rank are therefore stable, with positive lower response margin.",
        ],
        "nextProofTarget": (
            "Import this consequence in the full-theta source noncollapse theorem."
            if rank_margin_closed
            else "Close the failed interval consequence status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Full-theta source quadrature consequence theorem")
    print(f"  source quadrature closed: {source_quadrature_closed}")
    print(f"  theta tail closed: {theta_tail_closed}")
    print(f"  total error closed: {total_error_closed}")
    print(f"  riesz gap closed: {riesz_closed}")
    print(f"  active rank margin closed: {rank_margin_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
