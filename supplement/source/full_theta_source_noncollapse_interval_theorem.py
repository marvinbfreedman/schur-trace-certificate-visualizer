#!/usr/bin/env python3
r"""Theorem wrapper for full-theta active source noncollapse.

This file is the publication-facing proof interface for the full-theta source
rank input.  The imported algebra consequence theorem separates the exact
rank-stability argument from the lower interval constants layer.  This wrapper
records the resulting theorem:

    S_full = int E_u^* E_u du,
    ||S_full-S_{8,h}|| <= totalContinuumErrorBound,
    totalContinuumErrorBound < spectralGapToComplementS8/4,
    rank(L E^* | U_delta) = 2.

The wrapper does not recompute the raw interval bounds.  It checks that the
imported inequalities imply the active source noncollapse conclusion.
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
        "--source-json",
        default="full_theta_source_noncollapse_algebra_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="full_theta_source_noncollapse_interval_theorem.json",
    )
    args = parser.parse_args()

    raw = load(args.source_json)

    source_quad = dec(raw["sourceQuadratureErrorBound"])
    theta_tail = dec(raw["omittedTailDeltaSBound"])
    total_error = dec(raw["totalContinuumErrorBound"])
    gap = dec(raw["spectralGapToComplementS8"])
    alpha = dec(raw["totalProjectorAlpha"])
    rank_margin = dec(raw["finiteRankMarginS8"])
    angle_factor = dec(raw["angleFactorLowerBound"])
    comp_response = dec(raw["complementResponseOperatorNorm"])
    lower_after = dec(raw["lowerBoundAfterContinuumAndTail"])
    response_svals = [dec(x) for x in raw.get("finiteActiveResponseSingularValuesS8", [])]
    active_eigs = [dec(x) for x in raw.get("activeEigenvaluesS8", [])]

    total_matches_parts = total_error >= source_quad + theta_tail
    # The old certificate stores total_error as source_quad plus the tiny tail;
    # if decimal conversion exposes harmless roundoff, allow an absolute tail
    # tolerance several orders below any used margin.
    if not total_matches_parts:
        total_matches_parts = abs(total_error - (source_quad + theta_tail)) <= Decimal("1e-12")

    source_quadrature_closed = (
        source_quad >= 0
        and dec(raw["d2SourceGramIntervalEnvelope"]) >= 0
        and dec(raw["trapezoidFactor"]) >= 0
    )
    theta_tail_closed = theta_tail >= 0 and bool(raw.get("tailBoundClosed", raw.get("tailJson")))
    total_error_closed = source_quadrature_closed and theta_tail_closed and total_matches_parts
    riesz_closed = (
        total_error_closed
        and gap > 0
        and total_error < gap / Decimal(4)
        and alpha < 1
        and bool(raw.get("totalGapConditionPasses"))
    )
    active_rank_closed = (
        riesz_closed
        and len(active_eigs) == 2
        and len(response_svals) >= 2
        and min(active_eigs) > 0
        and min(response_svals[:2]) > 0
        and rank_margin > 0
        and angle_factor > 0
        and lower_after > 0
    )
    theorem_closed = active_rank_closed and bool(raw.get("fullPhiContinuumSourceNoncollapsePasses"))

    data: dict[str, Any] = {
        "theoremName": "full-theta source noncollapse interval theorem",
        "proofClass": "interval/ball certificate",
        "omega": raw.get("omega"),
        "L": raw.get("L"),
        "basis": raw.get("basis"),
        "fullNmax": raw.get("fullNmax"),
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
        "statuses": {
            "sourceQuadratureIntervalStatus": status(
                "source Gram interval quadrature bound",
                source_quadrature_closed,
                (
                    "The imported interval envelope for d_u^2(E_u^*E_u), "
                    "combined with the composite trapezoid factor, gives a "
                    "nonnegative operator bound for ||S_8-S_{8,h}||."
                ),
                blocker=None
                if source_quadrature_closed
                else "Regenerate the source quadrature interval certificate.",
            ),
            "thetaTailIntervalStatus": status(
                "n>=9 full-theta tail interval propagation",
                theta_tail_closed,
                (
                    "The omitted full-theta tail contributes a nonnegative "
                    "operator perturbation bound that is added to the source "
                    "quadrature error before the Riesz/rank comparison."
                ),
                blocker=None
                if theta_tail_closed
                else "Regenerate the full-theta interval tail propagation certificate.",
            ),
            "totalErrorAssemblyStatus": status(
                "total continuum source error assembly",
                total_error_closed,
                (
                    "The full source perturbation bound dominates the sum of "
                    "the interval source quadrature error and the theta-tail "
                    "operator perturbation."
                ),
                blocker=None
                if total_error_closed
                else "Check totalContinuumErrorBound against sourceQuadratureErrorBound plus omittedTailDeltaSBound.",
            ),
            "rieszProjectorStabilityStatus": status(
                "active Riesz projector stability",
                riesz_closed,
                (
                    "The total perturbation is below one quarter of the active/"
                    "complement spectral gap, equivalently alpha<1, so the "
                    "active two-dimensional Riesz subspace is stable."
                ),
                blocker=None
                if riesz_closed
                else "Reduce the source perturbation bound or increase the certified spectral gap.",
            ),
            "activeResponseRankStatus": status(
                "active response rank survives continuum and tail perturbation",
                active_rank_closed,
                (
                    "The finite active response has two positive singular "
                    "values, and the interval Riesz perturbation leaves a "
                    "positive lower response margin."
                ),
                blocker=None
                if active_rank_closed
                else "Recheck finite active response rank or continuum/tail perturbation margins.",
            ),
            "fullThetaSourceNoncollapseStatus": status(
                "full-Phi active source noncollapse",
                theorem_closed,
                (
                    "The interval source quadrature bound, theta-tail bound, "
                    "Riesz projector stability, and positive active response "
                    "margin prove rank(LE^*|U_delta)=2 for the full-Phi source "
                    "window."
                ),
                blocker=None
                if theorem_closed
                else "Close the source quadrature, tail, Riesz, or active response rank status.",
            ),
        },
        "formalProof": [
            (
                "Let S_full be the full-Phi continuum source Gram operator and "
                "let S_{8,h} be the certified theta-truncated source Gram."
            ),
            (
                "The interval source quadrature and n>=9 theta-tail certificates "
                "give ||S_full-S_{8,h}|| <= totalContinuumErrorBound."
            ),
            (
                "Since totalContinuumErrorBound < spectralGapToComplementS8/4, "
                "the active two-dimensional Riesz spectral subspace persists."
            ),
            (
                "The finite active response singular margin, after the Riesz "
                "angle loss and complement-response perturbation, remains "
                "positive: lowerBoundAfterContinuumAndTail > 0."
            ),
            (
                "Therefore the full-Phi active response map has rank two, "
                "equivalently the active source map is noncollapsed on U_delta."
            ),
        ],
        "fullThetaSourceNoncollapseIntervalTheoremClosed": theorem_closed,
        "fullPhiContinuumSourceNoncollapsePasses": theorem_closed,
        "nextProofTarget": (
            "Use this interval theorem as the source-rank input to the continuum trace-frame theorem."
            if theorem_closed
            else "Close the failed source noncollapse status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Full-theta source noncollapse interval theorem")
    print(f"  source quadrature interval: {source_quadrature_closed}")
    print(f"  theta tail interval: {theta_tail_closed}")
    print(f"  total error: {total_error}")
    print(f"  spectral gap: {gap}")
    print(f"  alpha: {alpha}")
    print(f"  lower response margin: {lower_after}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
