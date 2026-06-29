#!/usr/bin/env python3
r"""Continuum active trace range compatibility theorem ledger.

The endpoint Fredholm theorem reduces the adjoint Green endpoint condition to

    b(d_u) in Range(M).

The global Hilbert-space source of that compatibility is the active trace range
condition.  Let H be the A-Hilbert high block, let

    R f = (a -> Lambda_a(f)) on I,
    E_delta f = active source row of E f.

Then the exact continuum range criterion is

    E_delta in closure Range(R^*)      <=>      ker R subset ker E_delta.

This is the only statement strong enough to make the endpoint Fredholm
compatibility independent of sampled trace meshes.  Finite Galerkin injectivity
on the source-active subspace is useful evidence, but not a substitute for the
annihilator criterion: in the continuum proof one must rule out an
A-normalized closed-trace sequence with nonzero active source.

This ledger closes the functional-analytic equivalences and records the
remaining analytic unique-continuation/compactness input separately.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str | None) -> dict | None:
    if not path:
        return None
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


def summarize_active_range(active: dict | None) -> dict:
    if not active:
        return {}
    rows = active.get("rows", [])
    residuals = [
        float(row.get("rangeResidualRelative", 1.0))
        for row in rows
        if row.get("rangeResidualRelative") is not None
    ]
    ranks_match = all(
        int(row.get("traceRankOnActive", -1)) == int(row.get("activeDim", -2))
        for row in rows
    )
    coeff_ops = [
        float(row.get("coefficientOp", 0.0))
        for row in rows
        if row.get("coefficientOp") is not None
    ]
    return {
        "rows": len(rows),
        "ranksMatch": ranks_match,
        "maxRangeResidualRelative": max(residuals, default=None),
        "maxCoefficientOperator": max(coeff_ops, default=None),
        "source": "finite active block factorization E_active=C_sample R_global",
    }


def summarize_frame(frame: dict | None) -> dict:
    if not frame:
        return {}
    return {
        "rowCount": frame.get("rowCount"),
        "observedGammaFloor": frame.get("observedGammaFloor"),
        "allFrameMinsPositive": frame.get("allFrameMinsPositive"),
        "maxRangeResidualRelative": frame.get("maxRangeResidualRelative"),
        "maxDensityOperator": frame.get("maxDensityOperator"),
        "conditionalConclusion": frame.get("conditionalConclusion"),
    }


def summarize_high_block(high: dict | None) -> dict:
    if not high:
        return {}
    return {
        "tailEstimatePassesToContinuum": high.get("tailEstimatePassesToContinuum"),
        "conditionalTailEstimatePassage": high.get("conditionalTailEstimatePassage"),
        "ellipticTraceStatus": high.get("ellipticTraceStatus", {}),
        "moscoLimsupStatus": high.get("moscoLimsupStatus", {}),
        "moscoLiminfStatus": high.get("moscoLiminfStatus", {}),
        "sourceOperatorNormConvergenceStatus": high.get(
            "sourceOperatorNormConvergenceStatus", {}
        ),
    }


def summarize_unique_continuation(uc: dict | None) -> dict:
    if not uc:
        return {}
    return {
        "volterraSturmLagrangeIdentityStatus": uc.get(
            "volterraSturmLagrangeIdentityStatus", {}
        ),
        "distributionalInteriorGreenSourceStatus": uc.get(
            "distributionalInteriorGreenSourceStatus", {}
        ),
        "exactEndpointFredholmReductionStatus": uc.get(
            "exactEndpointFredholmReductionStatus", {}
        ),
        "formalGreenAnnihilationImplicationStatus": uc.get(
            "formalGreenAnnihilationImplicationStatus", {}
        ),
        "finiteOffBaseDerivativeRankUCStatus": uc.get(
            "finiteOffBaseDerivativeRankUCStatus", {}
        ),
        "sourceSideResponseNoncollapseStatus": uc.get(
            "sourceSideResponseNoncollapseStatus", {}
        ),
        "endpointCompatibilityForActiveRowsStatus": uc.get(
            "endpointCompatibilityForActiveRowsStatus", {}
        ),
        "closedTraceActiveUniqueContinuationStatus": uc.get(
            "closedTraceActiveUniqueContinuationStatus", {}
        ),
        "closedTraceActiveUniqueContinuationClosed": uc.get(
            "closedTraceActiveUniqueContinuationClosed"
        ),
        "remainingAnalyticGap": uc.get("remainingAnalyticGap"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-range-json", default="")
    parser.add_argument("--frame-json", default="")
    parser.add_argument("--quadrature-json", default="")
    parser.add_argument("--high-block-json", default="")
    parser.add_argument("--endpoint-range-json", default="adjoint_green_endpoint_range_interval_theorem.json")
    parser.add_argument(
        "--unique-continuation-json",
        default="closed_trace_active_unique_continuation_theorem.json",
    )
    parser.add_argument("--json-out", default="continuum_active_trace_range_theorem.json")
    args = parser.parse_args()

    active = load_optional(args.active_range_json)
    frame = load_optional(args.frame_json)
    quadrature = load_optional(args.quadrature_json)
    high = load_optional(args.high_block_json)
    endpoint = load_optional(args.endpoint_range_json)
    unique_continuation_json = load_optional(args.unique_continuation_json)

    active_summary = summarize_active_range(active)
    frame_summary = summarize_frame(frame)
    high_summary = summarize_high_block(high)
    unique_summary = summarize_unique_continuation(unique_continuation_json)

    annihilator = status(
        "Hilbert annihilator range criterion",
        True,
        (
            "For a densely defined closed trace map R between Hilbert spaces, "
            "closure Ran(R^*)=(ker R)^perp.  Hence E_active factors through "
            "the interval trace iff ker R is annihilated by E_active."
        ),
    )
    finite_factorization = status(
        "finite active block factorization",
        bool(
            active_summary
            and active_summary.get("ranksMatch")
            and (active_summary.get("maxRangeResidualRelative") or 1.0) < 1e-40
        ),
        (
            "Current Galerkin active blocks satisfy E_active=C_sample R_global "
            "to roundoff.  This is finite evidence for, not proof of, the "
            "continuum annihilator condition."
        ),
    )
    compactness_reduction = status(
        "compactness contradiction reduction",
        True,
        (
            "Failure of ker R subset ker E_active is equivalent to an "
            "A-normalized sequence with R f_n -> 0 and active source bounded "
            "below.  The high-block elliptic/Mosco theorem is the compactness "
            "mechanism needed to pass such a sequence to a continuum limit."
        ),
    )
    uc_imported = unique_summary.get("closedTraceActiveUniqueContinuationStatus", {})
    unique_continuation = status(
        "closed-trace active unique continuation",
        bool(uc_imported.get("closed")),
        (
            uc_imported.get("reason")
            if uc_imported
            else (
                "Still need the analytic theorem: a continuum limit with "
                "Lambda_a(f)=0 on I and nonzero active source cannot exist.  "
                "This must use the Volterra/Sturm Green structure, not finite "
                "sampling."
            )
        ),
    )
    continuum_range_closed = bool(unique_continuation["closed"] and annihilator["closed"])
    continuum_range = status(
        "continuum active trace range compatibility",
        continuum_range_closed,
        (
            "The closed-trace active unique-continuation theorem proves "
            "ker R subset ker E_active.  The Hilbert annihilator criterion "
            "then gives E_active in closure Ran(R^*)."
        ),
    )
    endpoint_compat = status(
        "endpoint Fredholm compatibility from active range",
        bool(continuum_range["closed"]),
        (
            "Since E_active lies in closure Ran(R^*), the endpoint vector "
            "b(d_u) annihilates ker(M^T), so b(d_u) lies in Range(M).  In "
            "the interval endpoint theorem this compatibility is also "
            "vacuous because the active endpoint map has full row rank."
        ),
    )

    data = {
        "theoremName": "continuum active trace range compatibility",
        "activeRangeJson": args.active_range_json if active else None,
        "frameJson": args.frame_json if frame else None,
        "quadratureJson": args.quadrature_json if quadrature else None,
        "highBlockJson": args.high_block_json if high else None,
        "endpointRangeJson": args.endpoint_range_json if endpoint else None,
        "uniqueContinuationJson": args.unique_continuation_json if unique_continuation_json else None,
        "hilbertAnnihilatorRangeCriterionStatus": annihilator,
        "finiteActiveBlockFactorizationStatus": finite_factorization,
        "compactnessContradictionReductionStatus": compactness_reduction,
        "closedTraceActiveUniqueContinuationStatus": unique_continuation,
        "continuumActiveTraceRangeCompatibilityStatus": continuum_range,
        "endpointFredholmCompatibilityFromActiveRangeStatus": endpoint_compat,
        "activeRangeEvidence": active_summary,
        "frameContinuumPassageEvidence": frame_summary,
        "quadratureStabilityEvidence": quadrature or {},
        "highBlockExhaustionTheorem": high_summary,
        "closedTraceActiveUniqueContinuationTheorem": unique_summary,
        "endpointRangeTheorem": endpoint or {},
        "exactEquivalences": [
            (
                "E_active in closure Ran(R^*) iff "
                "ker R subset ker E_active."
            ),
            (
                "Failure iff there are f_n with ||f_n||_A=1, "
                "||R f_n|| -> 0, and ||E_active f_n|| >= delta."
            ),
            (
                "If compactness gives f_n -> f locally and R f_n -> 0, then "
                "Lambda_a(f)=0 on I."
            ),
            (
                "The remaining unique-continuation lemma must show the same "
                "limit has E_active f=0, contradicting the lower active "
                "source bound."
            ),
        ],
        "formalProof": [
            (
                "Let R:H_A->L^2(I) be the closed interval trace operator and "
                "let E_active:H_A->U_delta be the active source map."
            ),
            (
                "By the closed-range/annihilator identity, a bounded active "
                "source row factors through R precisely when it vanishes on "
                "ker R."
            ),
            (
                "If vanishing fails, normalize a sequence in the bad set.  The "
                "commuted Sturm elliptic estimate plus Mosco high-block "
                "exhaustion is the compactness engine producing a nonzero "
                "closed-trace limit."
            ),
            (
                "The trace limit satisfies Lambda_a(f)=0 on the whole active "
                "interval.  The remaining analytic input is unique "
                "continuation in the Volterra/Sturm branch that forces the "
                "active source row of such a limit to vanish."
            ),
            (
                "That contradiction gives ker R subset ker E_active.  Hence "
                "E_active is in closure Ran(R^*), giving the continuum "
                "trace-to-source representation and the endpoint Fredholm "
                "compatibility b(d_u) in Range(M)."
            ),
        ],
        "remainingAnalyticGap": (
            None
            if continuum_range["closed"]
            else (
                unique_summary.get("remainingAnalyticGap")
                or (
                    "Prove the closed-trace active unique-continuation lemma under "
                    "the commuted Sturm/Mosco compactness hypotheses.  This is now "
                    "the specific continuum active-range bottleneck."
                )
            )
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum active trace range compatibility")
    print(f"  Hilbert annihilator criterion: {annihilator['closed']}")
    print(f"  finite active factorization: {finite_factorization['closed']}")
    if active_summary:
        print(
            "  finite max range residual: "
            f"{float(active_summary['maxRangeResidualRelative']):.12e}"
        )
    if frame_summary:
        print(
            "  observed frame gamma floor: "
            f"{float(frame_summary['observedGammaFloor']):.12e}"
        )
    print(f"  compactness reduction: {compactness_reduction['closed']}")
    print(f"  unique continuation closed: {unique_continuation['closed']}")
    print(f"  continuum active range closed: {continuum_range['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
