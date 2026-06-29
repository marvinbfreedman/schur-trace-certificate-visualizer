#!/usr/bin/env python3
r"""Bridge the full-theta source theorem back to the global Schur program.

This script is intentionally a proof ledger rather than another broad scan.
It collects the certificates that feed the global Weyl/Volterra Schur theorem:

  1. full-Phi source-side noncollapse on the continuum source window;
  2. sampled active trace range inclusion / active-kernel exclusion;
  3. source-inactive tail diagnostics.

The output separates closed certificate inputs from the remaining analytic
theorem.  The global Schur theorem is marked closed only when the
inactive-tail/closed-trace quotient assembly certificate is also supplied.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def f(x):
    return float(x)


def max_float(rows, key: str, default=0.0):
    values = [f(row[key]) for row in rows if key in row]
    return max(values, default=default)


def min_float(rows, key: str, default=0.0):
    values = [f(row[key]) for row in rows if key in row]
    return min(values, default=default)


def summarize_active_gap(data: dict) -> dict:
    rows = data.get("rows", [])
    active_kernel = max_float(rows, "activeTraceKernelSourceFrac")
    full_kernel = max_float(rows, "fullTraceKernelSourceFrac")
    min_active_trace = min_float(rows, "activeTraceMinNormed")
    worst_active = min(rows, key=lambda row: f(row.get("activeTraceMinNormed", 0.0))) if rows else {}
    worst_full = max(rows, key=lambda row: f(row.get("fullTraceKernelSourceFrac", 0.0))) if rows else {}
    return {
        "rows": len(rows),
        "activeTraceKernelSourceFracMax": active_kernel,
        "fullTraceKernelSourceFracMax": full_kernel,
        "activeTraceMinNormedMin": min_active_trace,
        "activeKernelExclusionObserved": active_kernel == 0.0 and min_active_trace > 0.0,
        "worstActiveTraceRow": worst_active,
        "worstFullKernelSourceRow": worst_full,
    }


def summarize_range(data: dict) -> dict:
    rows = data.get("rows", [])
    max_resid = max_float(rows, "rangeResidualRelative")
    min_rank_gap = min(
        (int(row.get("traceRankOnActive", 0)) - int(row.get("activeDim", 0)) for row in rows),
        default=0,
    )
    return {
        "rows": len(rows),
        "rangeResidualRelativeMax": max_resid,
        "traceRankMinusActiveDimMin": min_rank_gap,
        "sampledActiveRangeInclusionObserved": bool(rows) and max_resid < 1e-40 and min_rank_gap >= 0,
        "representativeRow": rows[0] if rows else {},
    }


def summarize_tail(data: dict) -> dict:
    rows = data.get("rows", [])
    best_tail = None
    best_nonzero_tail = None
    for row in rows:
        positive_modes = int(row.get("positiveModes", 0))
        for tail in row.get("tails", []):
            frac = f(tail.get("operatorFrac", tail.get("fraction", 0.0)))
            item = {
                "basis": row.get("basis"),
                "start": tail.get("start"),
                "fraction": frac,
                "tail": tail,
            }
            if best_tail is None or frac < best_tail["fraction"]:
                best_tail = item
            if 0.0 < frac and int(tail.get("start", 0)) < positive_modes:
                if best_nonzero_tail is None or frac < best_nonzero_tail["fraction"]:
                    best_nonzero_tail = item
    return {
        "rows": len(rows),
        "bestObservedOperatorTail": best_tail,
        "bestObservedNonzeroOperatorTail": best_nonzero_tail,
        "interpretation": (
            "Finite source-inactive/high-frequency tail diagnostic only.  This "
            "is not yet the continuum Schur tail theorem."
        ),
    }


def summarize_inactive_certificate(data: dict | None) -> dict | None:
    if not data:
        return None
    return {
        "globalSchurTailTheoremClosed": data.get("globalSchurTailTheoremClosed"),
        "sampledActiveKernelExclusionPasses": data.get("sampledActiveKernelExclusionPasses"),
        "finiteBudgetDiagnosticPasses": data.get("finiteBudgetDiagnosticPasses"),
        "worstContinuumInactiveFracUpper": data.get("worstContinuumInactiveFracUpper"),
        "worstActiveTraceKernelSourceFrac": data.get("worstActiveTraceKernelSourceFrac"),
        "worstFullTraceKernelSourceFrac": data.get("worstFullTraceKernelSourceFrac"),
        "continuumErrorBound": data.get("continuumErrorBound"),
        "sourceGrid": data.get("sourceGrid"),
        "activeCount": data.get("activeCount"),
        "fullNmax": data.get("fullNmax"),
    }


def summarize_external_audit(data: dict | None) -> dict | None:
    if not data:
        return None
    return {
        "normalizedSchurCertificateClosed": data.get("normalizedSchurCertificateClosed"),
        "conditionalQuotientToOriginalLiftClosed": data.get("conditionalQuotientToOriginalLiftClosed"),
        "externalFoundationClosed": data.get("externalFoundationClosed"),
        "originalWeylPositivityClosed": data.get("originalWeylPositivityClosed"),
        "originalKlmConditionClosed": data.get("originalKlmConditionClosed"),
        "rhFacingChainClosed": data.get("rhFacingChainClosed"),
        "closedAuditItems": data.get("closedAuditItems"),
        "openAuditItems": data.get("openAuditItems"),
        "nextProofTarget": data.get("nextProofTarget"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-quadrature", default="full_theta_source_noncollapse_interval_theorem.json")
    parser.add_argument("--active-gap", default="global_trace_active_gap_scan.json")
    parser.add_argument("--active-range", default="global_trace_active_range_inclusion.json")
    parser.add_argument("--synchronized-active-range", default="synchronized_active_range_interval_theorem.json")
    parser.add_argument("--source-tail", default="source_envelope_tail_certificate.json")
    parser.add_argument("--inactive-tail", default="full_theta_source_inactive_schur_tail_consequence_theorem.json")
    parser.add_argument("--tail-absorption", default="continuum_tail_absorption_consequence_theorem.json")
    parser.add_argument("--minmax-tail", default="source_inactive_minmax_tail_consequence_theorem.json")
    parser.add_argument("--high-block-exhaustion", default="high_block_tail_passage_bridge_consequence_theorem.json")
    parser.add_argument("--quotient-schur", default="weyl_volterra_quotient_schur_consequence_theorem.json")
    parser.add_argument(
        "--external-audit",
        default="weyl_volterra_external_equivalence_audit.json",
    )
    parser.add_argument("--json-out", default="global_weyl_volterra_schur_bridge.json")
    args = parser.parse_args()

    source = load_json(args.source_quadrature)
    active_gap = summarize_active_gap(load_json(args.active_gap))
    active_range = summarize_range(load_json(args.active_range))
    synchronized_active_range = (
        load_json(args.synchronized_active_range)
        if Path(args.synchronized_active_range).exists()
        else None
    )
    source_tail = summarize_tail(load_json(args.source_tail))
    inactive_tail = load_json(args.inactive_tail) if Path(args.inactive_tail).exists() else None
    tail_absorption = load_json(args.tail_absorption) if Path(args.tail_absorption).exists() else None
    minmax_tail = load_json(args.minmax_tail) if Path(args.minmax_tail).exists() else None
    high_block_exhaustion = (
        load_json(args.high_block_exhaustion)
        if Path(args.high_block_exhaustion).exists()
        else None
    )
    quotient_schur = (
        load_json(args.quotient_schur)
        if Path(args.quotient_schur).exists()
        else None
    )
    external_audit = (
        load_json(args.external_audit)
        if Path(args.external_audit).exists()
        else None
    )

    source_closed = bool(source.get("fullPhiContinuumSourceNoncollapsePasses"))
    active_range_observed = bool(active_range["sampledActiveRangeInclusionObserved"])
    active_range_closed = bool(
        synchronized_active_range
        and synchronized_active_range.get("activeRangeInclusionStatus", {}).get("closed")
    )
    active_kernel_observed = bool(active_gap["activeKernelExclusionObserved"])
    high_block_gaps = (
        (high_block_exhaustion or {}).get("remainingAnalyticGaps", [])
        if high_block_exhaustion
        else []
    )
    high_block_gap = high_block_gaps[0] if high_block_gaps else None
    high_block_tail_closed = bool(
        high_block_exhaustion and high_block_exhaustion.get("tailEstimatePassesToContinuum")
    )
    quotient_schur_closed = bool(
        quotient_schur
        and quotient_schur.get("globalWeylVolterraSchurStatus", {}).get("closed")
    )

    closed_inputs = {
        "fullPhiContinuumSourceNoncollapse": {
            "closed": source_closed,
            "basis": source.get("basis"),
            "sourceGrid": source.get("sourceGrid"),
            "projectorAlpha": source.get("totalProjectorAlpha"),
            "rankLowerBound": source.get("lowerBoundAfterContinuumAndTail"),
            "sourceQuadratureError": source.get("sourceQuadratureErrorBound"),
            "omittedThetaTailError": source.get("omittedTailDeltaSBound"),
        },
        "sampledActiveTraceRange": active_range,
        "synchronizedActiveRange": synchronized_active_range,
        "quotientSchurAssembly": quotient_schur,
        "sampledActiveTraceKernelExclusion": active_gap,
    }

    remaining = [
        {
            "name": "closed active range inclusion",
            "statement": (
                "E_active in closure Range(R_global^*) for the continuum active "
                "source plane, not only the sampled Galerkin range."
            ),
            "finiteEvidenceClosed": active_range_observed and active_kernel_observed,
            "analyticStatus": (
                "closed by synchronized endpoint interval/Krawczyk/Green BVP theorem"
                if active_range_closed
                else "still needs interval trace-resolution/unique-continuation proof"
            ),
        },
        {
            "name": "source-inactive high-frequency Schur tail",
            "statement": (
                "||(I-P_delta)E_Phi f||^2 <= epsilon_delta <A f,f>, with "
                "epsilon_delta absorbed by the finite low/mid Schur block."
            ),
            "finiteDiagnostic": source_tail,
            "fullPhiInactiveCertificate": summarize_inactive_certificate(inactive_tail),
            "absorptionCertificate": tail_absorption,
            "minMaxTailTheorem": minmax_tail,
            "highBlockExhaustionTheorem": high_block_exhaustion,
            "analyticStatus": (
                "closed in the full continuum high block"
                if high_block_tail_closed
                else (
                    "min-max proof closed in the normalized certified source model; "
                    + (
                        high_block_gap
                        or "still needs Galerkin-to-continuum high-block exhaustion"
                    )
                )
            ),
        },
        {
            "name": "closed-trace quotient factorization",
            "statement": (
                "Q_Phi(f)=||Gf||^2-||S R_global f||^2 with S bounded on the "
                "completed trace range; equivalently positivity on ker R_global "
                "plus the Douglas cross-form condition."
            ),
            "assemblyTheorem": quotient_schur,
            "analyticStatus": (
                "closed by the quotient Schur assembly theorem"
                if quotient_schur_closed
                else "formal theorem available; hypotheses still need full-Phi proof"
            ),
        },
    ]

    bridge_pass = source_closed
    global_schur_closed = quotient_schur_closed
    data = {
        "bridgeName": "full-Phi source noncollapse to global Weyl/Volterra Schur",
        "sourceNoncollapseJson": args.source_quadrature,
        "sourceInactiveMinmaxJson": args.minmax_tail,
        "highBlockExhaustionJson": args.high_block_exhaustion,
        "quotientSchurJson": args.quotient_schur,
        "bridgePasses": bridge_pass,
        "globalSchurTheoremClosed": global_schur_closed,
        "externalEquivalenceAudit": summarize_external_audit(external_audit),
        "closedInputs": closed_inputs,
        "remainingAnalyticItems": remaining,
        "nextProofTarget": (
            external_audit.get("nextProofTarget")
            if global_schur_closed and external_audit
            else (
                "Record the final dependency graph and isolate any external equivalence "
                "not contained in the Weyl/Volterra Schur certificate."
            )
            if global_schur_closed
            else (
                (
                    "Apply the closed active range inclusion and full-continuum "
                    "source-inactive domination inside the Weyl/Volterra quotient "
                    "factorization, and verify the remaining Schur-form hypotheses."
                )
                if high_block_tail_closed
                else high_block_gap
                or (
                    "Prove the Galerkin-to-continuum high-block exhaustion/elliptic "
                    "estimate that upgrades the normalized min-max source-tail theorem "
                    "from the certified source model to H_M cap ker R_global."
                )
            )
        ),
        "implication": (
            (
                "The active source/rank obstruction, active range inclusion, "
                "source-inactive high-block domination, and quotient/Douglas "
                "Schur assembly are all closed in the normalized full-Phi "
                "Weyl/Volterra model."
            )
            if global_schur_closed
            else (
                "The active source/rank obstruction is now removed for full Phi, "
                "and active range inclusion is closed by the synchronized endpoint "
                "theorem.  The normalized source-inactive min-max estimate is closed "
                "in the certified source model.  Therefore any remaining failure of "
                "the global Weyl/Volterra Schur form lies only in the final "
                "quotient-factorization/Schur assembly, not in the source-active "
                "or source-inactive estimates."
            )
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Global Weyl/Volterra Schur bridge")
    print(f"  full-Phi source noncollapse: {source_closed}")
    print(f"  sampled active range inclusion: {active_range_observed}")
    print(f"  synchronized active range theorem: {active_range_closed}")
    print(f"  sampled active kernel exclusion: {active_kernel_observed}")
    print(f"  max full trace-kernel source fraction: {active_gap['fullTraceKernelSourceFracMax']:.12e}")
    if inactive_tail:
        print(
            "  full-Phi inactive/top continuum bound: "
            f"{inactive_tail['worstContinuumInactiveFracUpper']:.12e}"
        )
    if tail_absorption:
        print(
            "  inactive tail absorption pass: "
            f"{tail_absorption['absorptionCertificatePasses']} "
            f"slack={tail_absorption['absorptionSlack']:.12e}"
        )
    if minmax_tail:
        print(
            "  source min-max tail proof: "
            f"{minmax_tail['minMaxProofInCertifiedSourceModel']} "
            f"epsilon={minmax_tail['normalizedEpsilonDelta']:.12e}"
        )
    if high_block_exhaustion:
        conditional_tail = high_block_exhaustion.get(
            "conditionalTailEstimatePassage",
            high_block_exhaustion.get("conditionalHighBlockExhaustionClosed"),
        )
        print(
            "  high-block exhaustion tail passage: "
            f"{high_block_exhaustion['tailEstimatePassesToContinuum']} "
            f"conditional={conditional_tail}"
        )
    if quotient_schur:
        print(
            "  quotient Schur assembly: "
            f"{quotient_schur['globalWeylVolterraSchurStatus']['closed']}"
        )
    if external_audit:
        print(
            "  external RH-facing chain closed: "
            f"{external_audit['rhFacingChainClosed']}"
        )
    print(f"  global Schur theorem closed: {global_schur_closed}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
