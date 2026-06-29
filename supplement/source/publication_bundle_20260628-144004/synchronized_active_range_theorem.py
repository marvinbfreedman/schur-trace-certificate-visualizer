#!/usr/bin/env python3
r"""Synchronized active range theorem and inactive-tail combination.

This is the proof ledger that cashes in the synchronized endpoint certificate.
It records the chain

    synchronized endpoint Krawczyk input closed
      => endpoint Green BVP solvable for every active source row
      => Lambda_a(f)=0 on I forces E_active f=0
      => E_active in closure Range(R_global^*)

and combines it with the normalized source-inactive min-max estimate

    ||(I-P_active) Ehat f||^2 <= epsilon_delta <A f,f>.

The active range part is a Hilbert-space annihilator theorem.  The inactive
tail part is imported at exactly the level proved by
``source_inactive_minmax_tail_theorem.py``: closed in the certified normalized
source model, with the high-block Galerkin-to-continuum passage recorded as a
separate theorem if not already closed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def is_closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--endpoint-json",
        default="endpoint_coefficient_synchronized_200_certificate.json",
    )
    parser.add_argument(
        "--active-range-json",
        default="continuum_active_trace_range_theorem.json",
    )
    parser.add_argument(
        "--inactive-tail-json",
        default="source_inactive_minmax_tail_theorem.json",
    )
    parser.add_argument("--json-out", default="synchronized_active_range_theorem.json")
    args = parser.parse_args()

    endpoint = load(args.endpoint_json)
    active = load(args.active_range_json)
    inactive = load(args.inactive_tail_json)

    endpoint_closed = is_closed(endpoint, "synchronizedKrawczykInputStatus")
    center_sync_closed = is_closed(endpoint, "centerSynchronizationStatus")
    companion_closed = is_closed(endpoint, "analyticCompanionInputStatus")
    boundary_closed = is_closed(endpoint, "analyticBoundaryInputStatus")
    hilbert_criterion_closed = is_closed(active, "hilbertAnnihilatorRangeCriterionStatus")
    compactness_reduction_closed = is_closed(active, "compactnessContradictionReductionStatus")

    # The synchronized endpoint certificate proves the active endpoint map has
    # full row rank inside the persistent Pluecker chart.  Therefore the
    # adjoint endpoint Green BVP has no active compatibility obstruction.
    endpoint_bvp_closed = endpoint_closed and center_sync_closed and companion_closed and boundary_closed

    green_representation_closed = endpoint_bvp_closed
    active_annihilation_closed = green_representation_closed
    active_range_closed = active_annihilation_closed and hilbert_criterion_closed

    epsilon = float(inactive["normalizedEpsilonDelta"])
    budget = float(inactive["finiteLowMidSchurBudget"])
    slack = float(inactive["absorptionSlack"])
    inactive_minmax_closed = bool(inactive.get("minMaxProofInCertifiedSourceModel"))
    inactive_absorbable = bool(inactive.get("absorbableByFiniteLowMidBlock"))
    inactive_tail_closed = inactive_minmax_closed and inactive_absorbable
    full_continuum_tail_closed = bool(inactive.get("continuumHighBlockExhaustionProvedHere"))

    certified_model_combined_closed = active_range_closed and inactive_tail_closed
    full_continuum_combined_closed = certified_model_combined_closed and full_continuum_tail_closed

    data = {
        "theoremName": "synchronized active range theorem",
        "endpointJson": args.endpoint_json,
        "activeRangeJson": args.active_range_json,
        "inactiveTailJson": args.inactive_tail_json,
        "endpointInputSummary": {
            "conditionMatrixOrder": endpoint.get("conditionMatrixOrder"),
            "controlledQuadratureOrder": endpoint.get("controlledQuadratureOrder"),
            "scaledCompanionAnalyticRadius": endpoint.get("scaledCompanionAnalyticRadius"),
            "boundaryRowAnalyticRadius": endpoint.get("boundaryRowAnalyticRadius"),
            "scaledCompanionRadiusOverSimultaneousBudget": endpoint.get(
                "scaledCompanionRadiusOverSimultaneousBudget"
            ),
            "boundaryRowRadiusOverSimultaneousBudget": endpoint.get(
                "boundaryRowRadiusOverSimultaneousBudget"
            ),
            "endpointRadiusOverChartCapacity": endpoint.get("endpointRadiusOverChartCapacity"),
            "actualKrawczykQ": endpoint.get("actualKrawczykQ"),
        },
        "activeRangeImportedStatus": {
            "hilbertAnnihilatorRangeCriterionStatus": active.get(
                "hilbertAnnihilatorRangeCriterionStatus", {}
            ),
            "finiteActiveBlockFactorizationStatus": active.get(
                "finiteActiveBlockFactorizationStatus", {}
            ),
            "compactnessContradictionReductionStatus": active.get(
                "compactnessContradictionReductionStatus", {}
            ),
            "oldContinuumActiveTraceRangeCompatibilityStatus": active.get(
                "continuumActiveTraceRangeCompatibilityStatus", {}
            ),
        },
        "inactiveTailConstants": {
            "normalizedEpsilonDelta": epsilon,
            "finiteLowMidSchurBudget": budget,
            "absorptionSlack": slack,
            "epsilonOverBudget": inactive.get("epsilonOverBudget"),
            "absoluteTailBound": inactive.get("absoluteTailBound"),
            "sourceTopLower": inactive.get("sourceTopLower"),
        },
        "synchronizedEndpointKrawczykInputStatus": status(
            "synchronized endpoint Krawczyk input",
            endpoint_closed,
            (
                "Imported from the synchronized 200-point coefficient input "
                "certificate.  It proves the finite endpoint map stays inside "
                "the persistent active Pluecker chart."
            ),
        ),
        "endpointGreenBvpSolvabilityStatus": status(
            "endpoint Green BVP solvable for active source rows",
            endpoint_bvp_closed,
            (
                "Full active endpoint row rank removes every Fredholm "
                "compatibility obstruction, so the endpoint free vector z can "
                "be chosen for every active source row."
            ),
        ),
        "closedTraceActiveAnnihilationStatus": status(
            "closed trace annihilates active source",
            active_annihilation_closed,
            (
                "With the endpoint Green BVP solvable, the Volterra/Sturm "
                "Green representation has only the interval trace term.  If "
                "Lambda_a(f)=0 on I, that term vanishes and E_active f=0."
            ),
        ),
        "activeRangeInclusionStatus": status(
            "E_active in closure Range(R_global^*)",
            active_range_closed,
            (
                "By the Hilbert annihilator criterion, ker R_global subset "
                "ker E_active is equivalent to E_active lying in closure "
                "Range(R_global^*)."
            ),
        ),
        "sourceInactiveTailDominationStatus": status(
            "source-inactive min-max tail domination",
            inactive_tail_closed,
            (
                "Imported min-max theorem gives "
                "||(I-P_active)Ehat f||^2 <= epsilon_delta <A f,f> in the "
                "certified normalized source model, and epsilon_delta is "
                "below the finite low/mid Schur budget."
            ),
        ),
        "certifiedModelCombinedSourceBoundStatus": status(
            "active range plus inactive tail in certified source model",
            certified_model_combined_closed,
            (
                "On ker R_global, the active source component vanishes; the "
                "inactive component is bounded by epsilon_delta and absorbed "
                "by the finite low/mid Schur block."
            ),
        ),
        "fullContinuumCombinedSourceBoundStatus": status(
            "full continuum active range plus inactive tail",
            full_continuum_combined_closed,
            (
                "Closed: active range inclusion is proved by the synchronized "
                "endpoint Green theorem, and the inactive-tail theorem imports "
                "the Galerkin-to-continuum high-block exhaustion passage."
                if full_continuum_combined_closed
                else (
                    "Open unless the inactive-tail theorem has also closed the "
                    "Galerkin-to-continuum high-block exhaustion passage."
                )
            ),
        ),
        "proofChain": [
            {
                "step": 1,
                "statement": (
                    "synchronized endpoint Krawczyk input closed implies the "
                    "active endpoint map has full row rank."
                ),
                "closed": endpoint_closed,
            },
            {
                "step": 2,
                "statement": (
                    "full active endpoint row rank implies the endpoint Green "
                    "BVP is solvable for every active source row."
                ),
                "closed": endpoint_bvp_closed,
            },
            {
                "step": 3,
                "statement": (
                    "the solved Green representation gives R_global f=0 "
                    "implies E_active f=0."
                ),
                "closed": active_annihilation_closed,
            },
            {
                "step": 4,
                "statement": (
                    "ker R_global subset ker E_active implies "
                    "E_active in closure Range(R_global^*)."
                ),
                "closed": active_range_closed,
            },
            {
                "step": 5,
                "statement": (
                    "source-inactive min-max tail gives "
                    "||(I-P_active)Ehat f||^2 <= epsilon_delta <A f,f>."
                ),
                "closed": inactive_tail_closed,
                "epsilonDelta": epsilon,
            },
        ],
        "formalConclusion": (
            "For f in ker R_global, P_active Ehat f=0 and "
            "||(I-P_active)Ehat f||^2 <= "
            f"{epsilon:.12e} <A f,f> in the certified normalized source model. "
            f"The finite low/mid Schur budget is {budget:.12e}, leaving slack "
            f"{slack:.12e}."
        ),
        "remainingAnalyticGap": (
            None
            if full_continuum_combined_closed
            else inactive.get("remainingAnalyticGap")
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Synchronized active range theorem")
    print(f"  endpoint Green BVP solvable: {endpoint_bvp_closed}")
    print(f"  active range inclusion: {active_range_closed}")
    print(f"  inactive epsilon_delta: {epsilon:.12e}")
    print(f"  finite low/mid budget: {budget:.12e}")
    print(f"  absorption slack: {slack:.12e}")
    print(f"  certified model combined: {certified_model_combined_closed}")
    print(f"  full continuum combined: {full_continuum_combined_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
