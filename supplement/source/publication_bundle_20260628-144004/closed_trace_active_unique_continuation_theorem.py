#!/usr/bin/env python3
r"""Closed-trace active unique-continuation theorem.

Target:

    Lambda_a(f)=0 on I  =>  E_active f=0

for closed-trace high-block limits.  This hardened ledger imports exact theorem
objects rather than sampled Green scans:

* symbolic Lagrange identity for P and P^*;
* symbolic distributional jump theorem for local source rows;
* endpoint interval/Krawczyk theorem giving full active row rank.

Together they construct, for every scalar active source row ell, a Green
coefficient K_ell with

    ell(f) = int_I K_ell(a)Lambda_a(f) da

on the completed trace domain.  Therefore closed trace implies active-source
annihilation.
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


def item_closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lagrange-json",
        default="trace_lagrange_adjoint_identity_consequence_theorem.json",
    )
    parser.add_argument(
        "--jump-json",
        default="adjoint_green_jump_conditions_consequence_theorem.json",
    )
    parser.add_argument(
        "--endpoint-range-json",
        default="adjoint_green_endpoint_range_interval_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="closed_trace_active_unique_continuation_theorem.json")
    args = parser.parse_args()

    lagrange = load(args.lagrange_json)
    jump = load(args.jump_json)
    endpoint = load(args.endpoint_range_json)

    lagrange_closed = bool(lagrange.get("theoremClosed")) and item_closed(
        lagrange, "lagrangeIdentityStatus"
    )
    jump_closed = bool(jump.get("theoremClosed")) and bool(
        jump.get("closedStatements", {}).get("interiorSourceSolved")
    )
    endpoint_reduction_closed = item_closed(
        endpoint, "exactFundamentalMatrixEndpointMapStatus"
    ) and item_closed(endpoint, "endpointFredholmRangeAlternativeStatus")
    endpoint_bvp_closed = item_closed(endpoint, "endpointGreenBvpSolvabilityStatus")
    endpoint_bound_closed = item_closed(
        endpoint, "uniformTraceDualBoundFromCompactnessStatus"
    )

    formal_green_closed = lagrange_closed and jump_closed and endpoint_reduction_closed
    theorem_closed = formal_green_closed and endpoint_bvp_closed and endpoint_bound_closed

    data = {
        "theoremName": "closed-trace active unique continuation",
        "proofClass": "analytic proof",
        "lagrangeJson": args.lagrange_json,
        "jumpJson": args.jump_json,
        "endpointRangeJson": args.endpoint_range_json,
        "volterraSturmLagrangeIdentityStatus": status(
            "Volterra/Sturm Lagrange identity",
            lagrange_closed,
            (
                "Imported from the symbolic Lagrange identity theorem: "
                "D_a B_P[K,f]=K P f-f P^*K."
            ),
        ),
        "distributionalInteriorGreenSourceStatus": status(
            "distributional interior Green source",
            jump_closed,
            (
                "Imported from the symbolic jump theorem: every local active "
                "source row has a unique jump vector, since the jump map is "
                "triangular with nonzero leading coefficient."
            ),
        ),
        "exactEndpointFredholmReductionStatus": status(
            "exact endpoint Fredholm reduction",
            endpoint_reduction_closed,
            (
                "The exact fundamental matrix gives the endpoint equation "
                "Mz=-b(ell), with the usual finite-dimensional Fredholm "
                "alternative."
            ),
        ),
        "endpointGreenBvpSolvabilityStatus": status(
            "endpoint Green BVP solvable for active rows",
            endpoint_bvp_closed,
            (
                "The endpoint interval/Krawczyk theorem proves full active "
                "row rank, so the endpoint compatibility condition is "
                "vacuous for active source rows."
            ),
        ),
        "formalGreenAnnihilationImplicationStatus": status(
            "formal Green annihilation implication",
            formal_green_closed and endpoint_bvp_closed,
            (
                "For each active scalar source row ell, the solved Green "
                "coefficient gives ell(f)=int_I K_ell(a)Lambda_a(f)da."
            ),
        ),
        "traceDualClosureStatus": status(
            "closed-domain trace-dual continuity",
            endpoint_bound_closed,
            (
                "The endpoint theorem supplies uniform trace-dual bounds for "
                "the selected Green family, so the identity passes from smooth "
                "test functions to the completed closed trace domain."
            ),
        ),
        "closedTraceActiveUniqueContinuationStatus": status(
            "closed-trace active unique continuation",
            theorem_closed,
            (
                "If Lambda_a(f)=0 on I, the Green representation annihilates "
                "every scalar active source row.  Therefore E_active f=0."
            ),
        ),
        "closedTraceActiveUniqueContinuationClosed": theorem_closed,
        "proof": [
            (
                "Let P f=Lambda_a(f).  The symbolic Lagrange identity gives "
                "D_a B_P[K,f]=K P f-f P^*K."
            ),
            (
                "For an active source row ell, the distributional jump theorem "
                "constructs a piecewise homogeneous K with P^*K=ell as a "
                "distribution."
            ),
            (
                "The endpoint interval theorem solves the endpoint concomitant "
                "equation for this K, because the active endpoint map has full "
                "row rank."
            ),
            (
                "Integrating Green's formula over I eliminates endpoint terms "
                "and gives ell(f)=int_I K(a)Lambda_a(f)da."
            ),
            (
                "Uniform trace-dual bounds pass this identity to the closed "
                "trace completion.  Hence R_global f=0 implies ell(f)=0 for "
                "each active row ell."
            ),
        ],
        "remainingAnalyticGap": None if theorem_closed else "One of the imported exact/interval theorem inputs is not closed.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Closed-trace active unique continuation")
    print(f"  Lagrange identity: {lagrange_closed}")
    print(f"  interior Green jumps: {jump_closed}")
    print(f"  endpoint Fredholm reduction: {endpoint_reduction_closed}")
    print(f"  endpoint Green BVP solvable: {endpoint_bvp_closed}")
    print(f"  trace-dual closure: {endpoint_bound_closed}")
    print(f"  closed-trace active UC closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
