#!/usr/bin/env python3
r"""Source-range Hardy/Green theorem ledger.

The target estimate is

    E_u^* E_u <= eta_u A,
    f in H_M cap ker R_global,

uniformly for u in the source window.  On the A-Hilbert high block this is
equivalent to a Green/Riesz representer theorem: every scalar component ell of
E_u has a representer g_ell with

    ell(f) = <g_ell,f>_A.

For the two-row source operator E_u, the optimal eta_u is the largest
eigenvalue of the 2x2 Gram matrix (<g_i,g_j>_A).

This ledger records that exact equivalence and imports the hardened continuum
Green-lift/active-range theorem.  The older finite Galerkin scans can still be
passed on the command line as diagnostics, but they are not proof inputs by
default.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.is_file():
        return None
    return json.loads(candidate.read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def max_float(rows: list[dict], key: str, default: float = 0.0) -> float:
    values = [float(row[key]) for row in rows if key in row and row[key] is not None]
    return max(values, default=default)


def summarize_green(green: dict | None) -> dict:
    if not green:
        return {}
    rows = green.get("rows", [])
    reps = green.get("representers", [])
    return {
        "basis": green.get("basis"),
        "constraints": green.get("constraints"),
        "cutoff": green.get("cutoff"),
        "highModes": green.get("highModes"),
        "maxRangeRelativeDefect": green.get("maxRangeRelativeDefect"),
        "maxEConstant": max_float(rows, "EConstant"),
        "maxDEConstant": max_float(rows, "dEConstant"),
        "maxCombinedConstant": max_float(rows, "combinedSumConstant"),
        "representerCount": len(reps),
        "rows": rows,
    }


def summarize_derivative(derivative: dict | None) -> dict:
    if not derivative:
        return {}
    return {
        "basis": derivative.get("basis"),
        "constraints": derivative.get("constraints"),
        "cutoff": derivative.get("cutoff"),
        "sourceWindow": [derivative.get("sourceMin"), derivative.get("sourceMax")],
        "sourceGrid": derivative.get("sourceGrid"),
        "mesh": derivative.get("mesh"),
        "meshRadius": derivative.get("meshRadius"),
        "maxSampleHighAbsorb": derivative.get("maxSampleHighAbsorb"),
        "minSampleFullAbsorb": derivative.get("minSampleFullAbsorb"),
        "maxSampleHighFullFrac": derivative.get("maxSampleHighFullFrac"),
        "maxAnalyticDerivativeHighAbsorb": derivative.get("maxAnalyticDerivativeHighAbsorb"),
        "maxAnalyticDerivativeHighFullFrac": derivative.get(
            "maxAnalyticDerivativeHighFullFrac"
        ),
        "analyticCoverHighAbsorb": derivative.get("analyticCoverHighAbsorb"),
        "analyticCoverHighFullFracVsMinFull": derivative.get(
            "analyticCoverHighFullFracVsMinFull"
        ),
        "worstDerivativeU": derivative.get("worstDerivativeU"),
        "finiteDifferenceChecks": derivative.get("finiteDifferenceChecks"),
    }


def summarize_aux(aux: dict | None) -> dict:
    if not aux:
        return {}
    return {
        "basis": aux.get("basis"),
        "constraints": aux.get("constraints"),
        "cutoff": aux.get("cutoff"),
        "sourceRangeAbsorb": aux.get("sourceRangeAbsorb"),
        "sourceRangeAbsorbFull": aux.get("sourceRangeAbsorbFull"),
        "sourceRangeAbsorbFrac": aux.get("sourceRangeAbsorbFrac"),
        "sourceRangeSize": aux.get("sourceRangeSize"),
        "perSource": aux.get("perSource"),
    }


def summarize_hardened(hardened: dict | None) -> dict:
    if not hardened:
        return {}
    statuses = hardened.get("statuses", {})
    return {
        "theoremName": hardened.get("theoremName"),
        "proofClass": hardened.get("proofClass"),
        "sourceWindow": hardened.get("sourceWindow"),
        "sourceRangeHardyGreenEstimateClosed": hardened.get(
            "sourceRangeHardyGreenEstimateClosed"
        ),
        "hardenedSourceRangeHardyGreenClosed": hardened.get(
            "hardenedSourceRangeHardyGreenClosed"
        ),
        "activeTraceControlStatus": statuses.get("activeTraceControlStatus"),
        "greenSourceControlStatus": statuses.get("greenSourceControlStatus"),
        "closedTraceHardyGreenEstimateStatus": statuses.get(
            "closedTraceHardyGreenEstimateStatus"
        ),
        "uniformRieszRepresenterStatus": statuses.get("uniformRieszRepresenterStatus"),
        "globalTraceSourceObservabilityStatus": statuses.get(
            "globalTraceSourceObservabilityStatus"
        ),
        "compactSourceOperatorStatus": statuses.get("compactSourceOperatorStatus"),
        "remainingAnalyticGap": hardened.get("remainingAnalyticGap"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hardened-json",
        default="source_range_hardy_green_hardened_theorem.json",
    )
    parser.add_argument("--green-json", default="")
    parser.add_argument("--derivative-json", default="")
    parser.add_argument("--aux-json", default="")
    parser.add_argument("--continuum-green-json", default="")
    parser.add_argument("--json-out", default="source_range_hardy_green_theorem.json")
    args = parser.parse_args()

    hardened = load_optional(args.hardened_json)
    green = load_optional(args.green_json)
    derivative = load_optional(args.derivative_json)
    aux = load_optional(args.aux_json)
    continuum_green = load_optional(args.continuum_green_json)
    hardened_summary = summarize_hardened(hardened)
    green_summary = summarize_green(green)
    derivative_summary = summarize_derivative(derivative)
    aux_summary = summarize_aux(aux)

    equivalence = status(
        "A-Riesz representer equivalence",
        True,
        (
            "In a Hilbert space, E_u^*E_u <= eta_u A is equivalent to the "
            "scalar source rows being A-bounded; the optimal eta_u is the top "
            "eigenvalue of the Gram matrix of their A-Riesz representers."
        ),
    )
    finite_representers = status(
        "finite Galerkin Green representers",
        bool(green_summary and float(green_summary.get("maxRangeRelativeDefect", 1.0)) < 1e-40),
        (
            "Finite closed-trace high-block representers solve the A-Riesz "
            "range equations to roundoff in the current Galerkin model."
        ),
    )
    source_window = status(
        "finite source-window derivative covering",
        bool(derivative_summary),
        (
            "The analytic u-derivative source rows give a finite-net covering "
            "bound over the source window; this is evidence for uniformity, "
            "not the continuum representer construction itself."
        ),
    )
    hardened_closed = bool(
        hardened and hardened.get("sourceRangeHardyGreenEstimateClosed")
    )
    continuum_green_closed = bool(
        continuum_green and continuum_green.get("continuumGreenRepresentersClosed")
    )
    continuum_closed = hardened_closed or continuum_green_closed
    if hardened_closed:
        continuum_reason = (
            "The hardened source-range theorem imports the completed "
            "Green-lift contraction and active-range interval theorem, giving "
            "A-bounded source rows and uniform A-Riesz representers without "
            "using scan-backed representer evidence."
        )
    elif continuum_green:
        continuum_reason = continuum_green.get("remainingAnalyticGap")
    else:
        continuum_reason = (
            "Need to construct g_{u,1}, g_{u,2} in the completed high-block "
            "A-space, prove ell_{u,k}(f)=<g_{u,k},f>_A for all closed-trace f, "
            "and bound the 2x2 representer Gram uniformly in u."
        )
    continuum_representers = status(
        "continuum closed-trace Green representers",
        continuum_closed,
        continuum_reason,
    )

    theorem_closed = continuum_representers["closed"]
    data = {
        "theoremName": "source-range Hardy/Green estimate",
        "hardenedJson": args.hardened_json if hardened else None,
        "greenJson": args.green_json if green else None,
        "derivativeJson": args.derivative_json if derivative else None,
        "auxiliaryJson": args.aux_json if aux else None,
        "continuumGreenJson": args.continuum_green_json if continuum_green else None,
        "rieszRepresenterEquivalenceStatus": equivalence,
        "finiteGreenRepresenterStatus": finite_representers,
        "sourceWindowDerivativeCoverStatus": source_window,
        "continuumGreenRepresenterStatus": continuum_representers,
        "hardenedSourceRangeTheorem": hardened_summary,
        "continuumGreenRepresenterTheorem": {
            "theoremName": continuum_green.get("theoremName"),
            "continuumGreenRepresentersClosed": continuum_green.get(
                "continuumGreenRepresentersClosed"
            ),
            "remainingAnalyticGap": continuum_green.get("remainingAnalyticGap"),
        }
        if continuum_green
        else None,
        "sourceRangeHardyGreenEstimateClosed": theorem_closed,
        "greenRepresenterCertificate": green_summary,
        "sourceWindowDerivativeCertificate": derivative_summary,
        "auxiliarySourceRangeCertificate": aux_summary,
        "formalProof": [
            (
                "Let H_hi be H_M cap ker R_global with inner product "
                "<f,g>_A=<Af,g>.  For fixed u, write E_u f=(ell_1(f),ell_2(f))."
            ),
            (
                "If ell_k(f)=<g_k,f>_A for g_k in H_hi, then "
                "||E_u f||^2 <= lambda_max((<g_i,g_j>_A)) ||f||_A^2 by "
                "Cauchy-Schwarz in the finite-dimensional representer span."
            ),
            (
                "Conversely, E_u^*E_u <= eta_u A makes each ell_k A-bounded, "
                "so the Riesz theorem supplies the representers.  Thus the "
                "estimate and Green representer theorem are equivalent."
            ),
            (
                "Uniform source-window control follows if u -> g_{u,k} is "
                "continuous in H_A and the representer Gram is bounded on "
                "the compact source window."
            ),
        ],
        "finiteEvidence": {
            "maxEConstant": green_summary.get("maxEConstant"),
            "maxDEConstant": green_summary.get("maxDEConstant"),
            "maxRangeRelativeDefect": green_summary.get("maxRangeRelativeDefect"),
            "analyticCoverHighAbsorb": derivative_summary.get("analyticCoverHighAbsorb"),
            "analyticCoverHighFullFracVsMinFull": derivative_summary.get(
                "analyticCoverHighFullFracVsMinFull"
            ),
            "sourceRangeAbsorbFrac": aux_summary.get("sourceRangeAbsorbFrac"),
        },
        "remainingAnalyticGap": None if theorem_closed else continuum_representers["reason"],
        "interpretation": (
            "The source-range Hardy/Green theorem is closed by the hardened "
            "continuum Green-lift/active-range theorem.  Finite representer "
            "and source-window scans remain useful diagnostics, but are not "
            "proof dependencies of this ledger."
            if theorem_closed
            else (
                "The source-range Hardy/Green theorem is reduced to a precise "
                "continuum Green representer construction.  Optional finite "
                "diagnostics do not close the theorem without the hardened "
                "continuum input."
            )
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Source-range Hardy/Green theorem")
    print(f"  Riesz equivalence closed: {equivalence['closed']}")
    print(f"  hardened continuum theorem closed: {hardened_closed}")
    print(f"  finite representers closed: {finite_representers['closed']}")
    if green_summary:
        print(f"  max finite E constant: {green_summary['maxEConstant']:.12e}")
        print(f"  max range defect: {float(green_summary['maxRangeRelativeDefect']):.12e}")
    if derivative_summary:
        print(
            "  analytic cover high/full fraction: "
            f"{float(derivative_summary['analyticCoverHighFullFracVsMinFull']):.12e}"
        )
    print(f"  continuum Green representers closed: {continuum_representers['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
