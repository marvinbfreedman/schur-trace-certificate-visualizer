#!/usr/bin/env python3
r"""Symbolic algebra consequence for full-theta source noncollapse.

The interval layer supplies certified constants for the source Gram
perturbation, theta-tail perturbation, active/complement spectral gap, and
finite active response margin.  This theorem contains only the deterministic
algebra:

    perturbation < gap/4
    positive finite response margin
        => stable two-dimensional active Riesz subspace
        => rank(LE^*|U_delta)=2 for the full-Phi source model.

Keeping this implication separate lets the audit distinguish the machine
interval enclosure from the exact rank-stability argument.
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
        "--constants-json",
        default="full_theta_source_quadrature_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="full_theta_source_noncollapse_algebra_consequence_theorem.json",
    )
    args = parser.parse_args()

    constants = load(args.constants_json)

    total_error = dec(constants["totalContinuumErrorBound"])
    gap = dec(constants["spectralGapToComplementS8"])
    alpha = dec(constants["totalProjectorAlpha"])
    lower_after = dec(constants["lowerBoundAfterContinuumAndTail"])
    rank_margin = dec(constants["finiteRankMarginS8"])
    active_eigs = [dec(x) for x in constants.get("activeEigenvaluesS8", [])]
    response_svals = [
        dec(x) for x in constants.get("finiteActiveResponseSingularValuesS8", [])
    ]

    constants_closed = bool(constants.get("sourceQuadratureConsequenceClosed"))
    gap_condition = gap > 0 and total_error < gap / Decimal(4) and alpha < 1
    finite_active_rank = (
        len(active_eigs) == 2
        and len(response_svals) >= 2
        and min(active_eigs) > 0
        and min(response_svals[:2]) > 0
        and rank_margin > 0
    )
    noncollapse_closed = constants_closed and gap_condition and finite_active_rank and lower_after > 0

    data: dict[str, Any] = {
        "theoremName": "full-theta source noncollapse algebra consequence theorem",
        "proofClass": "symbolic identity",
        "activeDimension": len(active_eigs),
        "activeEigenvaluesS8": constants.get("activeEigenvaluesS8", []),
        "finiteActiveResponseSingularValuesS8": constants.get(
            "finiteActiveResponseSingularValuesS8", []
        ),
        "sourceQuadratureErrorBound": constants.get("sourceQuadratureErrorBound"),
        "omittedTailDeltaSBound": constants.get("omittedTailDeltaSBound"),
        "totalContinuumErrorBound": constants.get("totalContinuumErrorBound"),
        "spectralGapToComplementS8": constants.get("spectralGapToComplementS8"),
        "totalProjectorAlpha": constants.get("totalProjectorAlpha"),
        "finiteRankMarginS8": constants.get("finiteRankMarginS8"),
        "angleFactorLowerBound": constants.get("angleFactorLowerBound"),
        "complementResponseOperatorNorm": constants.get("complementResponseOperatorNorm"),
        "lowerBoundAfterContinuumAndTail": constants.get("lowerBoundAfterContinuumAndTail"),
        "d2SourceGramIntervalEnvelope": constants.get("d2SourceGramIntervalEnvelope"),
        "trapezoidFactor": constants.get("trapezoidFactor"),
        "tailBoundClosed": bool(constants.get("tailBoundClosed")),
        "totalGapConditionPasses": bool(constants.get("totalGapConditionPasses")),
        "sourceQuadratureConsequenceClosed": constants_closed,
        "statuses": {
            "constantsInputStatus": status(
                "closed source quadrature constants",
                constants_closed,
                "The imported constants theorem has closed the source perturbation and theta-tail bounds.",
                blocker=None if constants_closed else "Close the source quadrature constants consequence.",
            ),
            "gapStabilityStatus": status(
                "rank stability gap inequality",
                gap_condition,
                "The perturbation is less than one quarter of the active/complement gap, so the active Riesz subspace is stable.",
                blocker=None if gap_condition else "Improve the perturbation/gap margin.",
            ),
            "finiteActiveRankStatus": status(
                "finite active response rank",
                finite_active_rank,
                "The two finite active response singular values and rank margin are positive.",
                blocker=None if finite_active_rank else "Certify the finite active response rank margin.",
            ),
            "fullThetaSourceNoncollapseStatus": status(
                "full-Phi source noncollapse algebra consequence",
                noncollapse_closed,
                "The gap-stability inequality and positive response margin imply rank(LE^*|U_delta)=2.",
                blocker=None if noncollapse_closed else "Close the constants, gap, or finite-rank status.",
            ),
        },
        "formalProof": [
            "The source perturbation bound is below one quarter of the active/complement gap.",
            "Davis-Kahan/Riesz stability keeps a two-dimensional active source subspace.",
            "The finite active response has positive smallest singular value and positive rank margin.",
            "The transported active response therefore remains rank two for the full-Phi source model.",
        ],
        "fullThetaSourceNoncollapseAlgebraClosed": noncollapse_closed,
        "fullPhiContinuumSourceNoncollapsePasses": noncollapse_closed,
        "nextProofTarget": (
            "Import this symbolic implication in the full-theta source noncollapse theorem."
            if noncollapse_closed
            else "Close the failed algebra consequence status above."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Full-theta source noncollapse algebra consequence theorem")
    print(f"  constants closed: {constants_closed}")
    print(f"  gap condition: {gap_condition}")
    print(f"  finite active rank: {finite_active_rank}")
    print(f"  theorem closed: {noncollapse_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
