#!/usr/bin/env python3
r"""Commuted Sturm / source-range elliptic trace theorem ledger.

Earlier notes proposed proving a local principal-coefficient identity

    S_10(f)=int a_10(s)|D^10 f|^2 + lower terms + trace squares,
    a_10 > 0.

That target is not valid for the endpoint compact kernel itself: high-frequency
interior packets make the compact commuted-kernel form vanish relative to the
Sobolev graph norm.  The corrected theorem is a source-range Hardy/Green
estimate on the closed-trace high block:

    E_u^* E_u <= eta_u A,        f in H_M cap ker R_global,

uniformly over the source window, with the source-range residual small enough
to be absorbed by the already certified finite low/mid Schur block.

This ledger records the finite evidence, the compact-kernel obstruction, and
the exact analytic sublemma still required for the high-block exhaustion
theorem.
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


def min_subunit_row(dominance: dict | None) -> dict | None:
    if not dominance:
        return None
    rows = [row for row in dominance.get("rows", []) if row.get("subunit")]
    if not rows:
        return None
    return min(rows, key=lambda row: float(row["product"]))


def compact_obstruction_summary(obstruction: dict | None) -> dict:
    if not obstruction:
        return {}
    rows = obstruction.get("rows", [])
    ratios = [float(row["ratio"]) for row in rows if "ratio" in row]
    return {
        "order": obstruction.get("order"),
        "support": obstruction.get("support"),
        "frequencyCount": len(rows),
        "maxRatio": max(ratios) if ratios else None,
        "minRatio": min(ratios) if ratios else None,
        "interpretation": obstruction.get("interpretation"),
    }


def source_range_summary(aux: dict | None) -> dict:
    if not aux:
        return {}
    per_source = aux.get("perSource", [])
    fracs = [float(row["sourceAbsorbFrac"]) for row in per_source if "sourceAbsorbFrac" in row]
    etas = [float(row["sourceAbsorb"]) for row in per_source if "sourceAbsorb" in row]
    return {
        "basis": aux.get("basis"),
        "constraints": aux.get("constraints"),
        "cutoff": aux.get("cutoff"),
        "positiveModes": aux.get("positiveModes"),
        "sobolevRouteProduct": aux.get("sobolevRouteProduct"),
        "sourceRangeAbsorb": aux.get("sourceRangeAbsorb"),
        "sourceRangeAbsorbFull": aux.get("sourceRangeAbsorbFull"),
        "sourceRangeAbsorbFrac": aux.get("sourceRangeAbsorbFrac"),
        "maxPerSourceAbsorb": max(etas) if etas else None,
        "maxPerSourceFrac": max(fracs) if fracs else None,
        "perSource": per_source,
        "interpretation": aux.get("interpretation"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dominance-json", default="lagrange_commuted_dominance_summary.json")
    parser.add_argument("--obstruction-json", default="commuted_compact_obstruction.json")
    parser.add_argument("--aux-json", default="aux_regularizer_certificate.json")
    parser.add_argument("--source-range-json", default="source_range_hardy_green_theorem.json")
    parser.add_argument(
        "--compact-proof-json",
        default="high_block_compact_exhaustion_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="commuted_sturm_elliptic_trace_theorem.json")
    args = parser.parse_args()

    dominance = load_optional(args.dominance_json)
    obstruction = load_optional(args.obstruction_json)
    aux = load_optional(args.aux_json)
    source_range = load_optional(args.source_range_json)
    compact_proof = load_optional(args.compact_proof_json)
    best_subunit = min_subunit_row(dominance)
    obstruction_data = compact_obstruction_summary(obstruction)
    source_data = source_range_summary(aux)

    compact_kernel_ellipticity = status(
        "compact commuted-kernel Sobolev ellipticity",
        False,
        (
            "False as a continuum theorem.  Compactly supported oscillatory "
            "packets in ker R make the compact commuted-kernel form vanish "
            "relative to W^{m,2}; the notes and obstruction certificate rule "
            "out proving coercivity from this object alone."
        ),
    )
    finite_commuted_domination = status(
        "finite commuted domination",
        bool(best_subunit),
        (
            "Finite Galerkin evidence shows subunit domination for suitable "
            "cutoff/order pairs, but this uses the compact commuted matrix and "
            "does not itself prove the continuum trace theorem."
        ),
    )
    source_range_closed = bool(
        source_range and source_range.get("sourceRangeHardyGreenEstimateClosed")
    )
    source_range_reason = (
        source_range.get("remainingAnalyticGap")
        if source_range
        else (
            "Need a direct continuum proof of E_u^*E_u <= eta_u A on "
            "H_M cap ker R_global, uniformly over the source window, using the "
            "closed-trace Green identity rather than full Sobolev coercivity."
        )
    )
    source_range_hardy_green = status(
        "continuum source-range Hardy/Green estimate",
        source_range_closed,
        source_range_reason,
    )
    compact_source_exhaustion = status(
        "compact-source high-block exhaustion implication",
        bool(compact_proof and compact_proof.get("conditionalHighBlockExhaustionClosed")),
        (
            "The corrected high-block proof is a compact-source/Mosco theorem, "
            "not a raw compact-kernel Sobolev coercivity theorem.  It closes the "
            "functional-analytic implication conditional on the continuum "
            "trace-frame lower bound."
        ),
    )
    endpoint_trace_boundedness = status(
        "endpoint source-row trace boundedness",
        bool(compact_proof and compact_proof.get("tailEstimatePassesToContinuum")),
        (
            "Follows once the compact-source exhaustion theorem is unconditional; "
            "the current remaining input is the continuum trace-frame lower "
            "bound for sampled trace correction."
        ),
    )

    theorem_closed = (
        source_range_hardy_green["closed"] and endpoint_trace_boundedness["closed"]
    ) or bool(compact_proof and compact_proof.get("tailEstimatePassesToContinuum"))
    data = {
        "theoremName": "commuted Sturm/source-range elliptic trace theorem",
        "dominanceJson": args.dominance_json if dominance else None,
        "obstructionJson": args.obstruction_json if obstruction else None,
        "auxiliaryJson": args.aux_json if aux else None,
        "sourceRangeJson": args.source_range_json if source_range else None,
        "compactProofJson": args.compact_proof_json if compact_proof else None,
        "compactKernelEllipticityStatus": compact_kernel_ellipticity,
        "finiteCommutedDominationStatus": finite_commuted_domination,
        "sourceRangeHardyGreenStatus": source_range_hardy_green,
        "sourceRangeHardyGreenTheorem": source_range,
        "compactSourceExhaustionStatus": compact_source_exhaustion,
        "compactSourceExhaustionProofSummary": (
            {
                "theoremName": compact_proof.get("theoremName"),
                "proofClass": compact_proof.get("proofClass"),
                "conditionalHighBlockExhaustionClosed": bool(
                    compact_proof.get("conditionalHighBlockExhaustionClosed")
                ),
                "tailEstimatePassesToContinuum": bool(
                    compact_proof.get("tailEstimatePassesToContinuum")
                ),
            }
            if compact_proof
            else None
        ),
        "endpointTraceBoundednessStatus": endpoint_trace_boundedness,
        "ellipticTraceEstimateClosed": theorem_closed,
        "highBlockExhaustionInputClosed": theorem_closed,
        "bestFiniteSubunitDominance": best_subunit,
        "compactObstruction": obstruction_data,
        "sourceRangeCertificate": source_data,
        "correctedTheorem": (
            "For the closed continuum trace space and after removing the finite "
            "low/mid A-modes, prove E_u^*E_u <= eta_u A uniformly in u on the "
            "source window.  This source-range Hardy/Green estimate is the "
            "elliptic trace input needed by the high-block exhaustion theorem."
        ),
        "invalidRoute": (
            "Do not try to prove raw local a_10(s)>0 Sobolev ellipticity for "
            "the compact endpoint commuted-kernel form; the compact packet "
            "obstruction disproves that route."
        ),
        "formalConsequencesIfClosed": [
            "Endpoint source rows are A-bounded on H_M cap ker R_global.",
            "The normalized source operator Ehat_Phi is compact/finite-rank in the A-Hilbert source model.",
            "P_N Ehat^*Ehat P_N -> Ehat^*Ehat follows once A-graph Galerkin/Mosco convergence is also proved.",
            "The high-block exhaustion ledger can then close sourceOperatorNormConvergenceStatus after Mosco is supplied.",
        ],
        "remainingAnalyticGap": (
            None
            if theorem_closed
            else (
                (compact_proof or {}).get("remainingAnalyticGap")
                or (
                    "Prove the continuum source-range Hardy/Green estimate directly "
                    "from the closed-trace Green identity and source-window coefficient "
                    "bounds."
                )
            )
        ),
        "interpretation": (
            "The corrected elliptic trace target has been sharpened: the "
            "compact-kernel Sobolev route is invalid, while the compact-source "
            "exhaustion implication is proved conditionally.  The remaining "
            "non-formal input is the continuum trace-frame lower bound."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Commuted Sturm/source-range elliptic trace theorem")
    print(f"  compact-kernel Sobolev route closed: {compact_kernel_ellipticity['closed']}")
    print(f"  finite commuted domination observed: {finite_commuted_domination['closed']}")
    if best_subunit:
        print(
            "  best finite subunit: "
            f"M={best_subunit['cutoff']} m={best_subunit['order']} "
            f"product={best_subunit['product']:.12e}"
        )
    if source_data:
        print(
            "  source-range high/full fraction: "
            f"{float(source_data['sourceRangeAbsorbFrac']):.12e}"
        )
    print(f"  continuum source-range Hardy/Green closed: {source_range_hardy_green['closed']}")
    print(f"  compact-source exhaustion implication: {compact_source_exhaustion['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
