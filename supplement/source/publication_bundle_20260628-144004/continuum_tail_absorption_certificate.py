#!/usr/bin/env python3
r"""Continuum tail absorption certificate.

This is the last algebraic step after the full-Phi source-inactive certificate.
It does not pretend to prove the continuum Hardy/Schur tail theorem.  Instead
it records the exact conditional theorem:

    if ||(I-P_delta)E_Phi f||^2 <= epsilon_delta <A f,f>
       on H_M cap ker R_global,

and if epsilon_delta is bounded by the certified full-Phi inactive constant,
then the finite low/mid Schur block absorbs the tail whenever

    epsilon_delta < beta_low_mid.

The script loads the full-Phi inactive-tail certificate and the finite Schur
tail-budget diagnostic, checks the numerical absorption inequality, and emits
the theorem ledger.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inactive-json",
        default="full_theta_source_inactive_schur_tail_certificate.json",
    )
    parser.add_argument("--json-out", default="continuum_tail_absorption_certificate.json")
    args = parser.parse_args()

    inactive = load(args.inactive_json)
    epsilon = float(inactive["worstContinuumInactiveFracUpper"])
    budget = float(inactive["finiteSchurTailBudgetDiagnostic"]["fraction"])
    margin = budget - epsilon
    ratio = epsilon / budget if budget else float("inf")
    absorption_pass = margin > 0

    row = inactive["rows"][0] if inactive.get("rows") else {}
    data = {
        "theoremName": "source-inactive continuum tail absorption",
        "inactiveJson": args.inactive_json,
        "continuumTailEstimateProved": False,
        "absorptionCertificatePasses": absorption_pass,
        "epsilonDeltaUpper": epsilon,
        "finiteLowMidSchurBudget": budget,
        "absorptionSlack": margin,
        "epsilonOverBudget": ratio,
        "sourceGrid": inactive.get("sourceGrid"),
        "basis": row.get("basis"),
        "globalTraceRatio": row.get("ratio"),
        "activeTraceKernelSourceFrac": inactive.get("worstActiveTraceKernelSourceFrac"),
        "fullTraceKernelSourceFrac": inactive.get("worstFullTraceKernelSourceFrac"),
        "inactiveTraceKernelSourceFrac": row.get("inactiveTraceKernelSourceFrac"),
        "continuumInactiveTopUpper": inactive.get("worstContinuumInactiveFracUpper"),
        "finiteInactiveTop": row.get("finiteInactiveFracOfTop"),
        "continuumErrorBound": inactive.get("continuumErrorBound"),
        "conditionalTheorem": (
            "If the continuum high-frequency Hardy/Schur estimate "
            "||(I-P_delta)E_Phi f||^2 <= epsilon_delta <A f,f> is proved on "
            "H_M cap ker R_global with epsilon_delta <= epsilonDeltaUpper, "
            "then the finite low/mid Schur block absorbs this source-inactive "
            "tail because epsilonDeltaUpper < finiteLowMidSchurBudget."
        ),
        "remainingAnalyticGap": (
            "Prove the continuum Hardy/Schur estimate itself, not just the "
            "finite spectral inactive bound.  The numerical absorption margin "
            "is positive once that theorem supplies the stated epsilon_delta."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuum tail absorption certificate")
    print(f"  epsilon_delta upper: {epsilon:.12e}")
    print(f"  finite low/mid Schur budget: {budget:.12e}")
    print(f"  slack: {margin:.12e}")
    print(f"  epsilon/budget: {ratio:.12e}")
    print(f"  absorption pass: {absorption_pass}")
    print("  continuum tail estimate proved: False")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
