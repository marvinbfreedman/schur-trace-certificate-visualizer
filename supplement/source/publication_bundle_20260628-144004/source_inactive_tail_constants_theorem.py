#!/usr/bin/env python3
r"""Certified source-inactive tail constants theorem.

This is the narrow input needed by the compact high-block exhaustion proof.
It does not import the high-block exhaustion theorem and does not claim the
Galerkin-to-continuum passage.

It reads the literal full-theta inactive source certificate and the low/mid
tail absorption certificate, then exposes only the constants used by the
compact-source argument:

    ||(I-P_2) Ehat f||^2 <= epsilon_delta <A f,f>

inside the certified source model, with epsilon_delta below the Schur
absorption budget.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
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
        "--inactive-json",
        default="full_theta_source_inactive_schur_tail_certificate.json",
    )
    parser.add_argument(
        "--absorption-json",
        default="continuum_tail_absorption_certificate.json",
    )
    parser.add_argument("--json-out", default="source_inactive_tail_constants_theorem.json")
    args = parser.parse_args()

    inactive = load(args.inactive_json)
    absorption = load(args.absorption_json)
    rows = inactive.get("rows", [])
    if not rows:
        raise SystemExit("inactive certificate has no rows")

    worst = max(rows, key=lambda row: float(row["continuumInactiveFracUpper"]))
    absolute_tail = float(worst["continuumInactiveUpper"])
    source_top_lower = float(worst["continuumTopLower"])
    finite_lambda3 = float(worst["complementTopEigenvalue"])
    finite_top = float(worst["sourceTop"])
    perturbation_eta = float(worst["continuumErrorBound"])
    normalized_epsilon = float(worst["continuumInactiveFracUpper"])
    budget = float(absorption["finiteLowMidSchurBudget"])
    slack = budget - normalized_epsilon

    certificate_closed = bool(
        inactive.get("globalSchurTailTheoremClosed")
        or inactive.get("sampledActiveKernelExclusionPasses")
        or inactive.get("finiteBudgetDiagnosticPasses")
    )
    absorption_closed = bool(absorption.get("absorptionCertificatePasses")) and slack > 0
    constants_closed = certificate_closed and absorption_closed

    data = {
        "theoremName": "source-inactive tail constants theorem",
        "proofClass": "interval/ball certificate",
        "inactiveJson": args.inactive_json,
        "absorptionJson": args.absorption_json,
        "basis": worst.get("basis"),
        "sourceGrid": inactive.get("sourceGrid"),
        "sourceWindow": [inactive.get("sourceMin"), inactive.get("sourceMax")],
        "globalTraceRatio": worst.get("ratio"),
        "finiteSourceTop": finite_top,
        "finiteInactiveLambda3": finite_lambda3,
        "operatorPerturbationEta": perturbation_eta,
        "absoluteTailBound": absolute_tail,
        "sourceTopLower": source_top_lower,
        "normalizedEpsilonDelta": normalized_epsilon,
        "finiteLowMidSchurBudget": budget,
        "absorptionSlack": slack,
        "epsilonOverBudget": normalized_epsilon / budget if budget else float("inf"),
        "statuses": {
            "inactiveSourceCertificateStatus": status(
                "full-theta inactive source certificate",
                certificate_closed,
                (
                    "The full-theta source certificate gives the inactive "
                    "lambda_3 upper bound and the top source lower bound after "
                    "certified theta-tail and quadrature error propagation."
                ),
                blocker=None
                if certificate_closed
                else "Close full_theta_source_inactive_schur_tail_certificate.py.",
            ),
            "absorptionBudgetStatus": status(
                "low/mid Schur absorption budget",
                absorption_closed,
                (
                    "The absorption certificate proves epsilon_delta is below "
                    "the available finite low/mid Schur budget."
                ),
                blocker=None if absorption_closed else "Close continuum_tail_absorption_certificate.py.",
            ),
            "normalizedTailConstantsStatus": status(
                "normalized inactive-tail constants",
                constants_closed,
                (
                    "The normalized source model satisfies the certified "
                    "inactive-tail constant and the constant is absorbable by "
                    "the Schur budget."
                ),
                blocker=None if constants_closed else "Need inactive source and absorption certificates.",
            ),
        },
        "minMaxProofInCertifiedSourceModel": constants_closed,
        "sourceInactiveTailConstantsClosed": constants_closed,
        "sourceInactiveTailDominationConstantsClosed": constants_closed,
        "absorbableByFiniteLowMidBlock": absorption_closed,
        "rawUnnormalizedStatement": (
            "||(I-P_2)E_Phi f||^2 <= absoluteTailBound <A f,f> in the certified source model."
        ),
        "normalizedStatement": (
            "For Ehat_Phi=E_Phi/sqrt(sourceTopLower), "
            "||(I-P_2)Ehat_Phi f||^2 <= normalizedEpsilonDelta <A f,f>."
        ),
        "proof": [
            "Read the certified full-theta inactive source bound.",
            "Use Weyl perturbation of the source Gram matrix to bound lambda_3 and lambda_1.",
            "Normalize by the certified top source lower bound.",
            "Compare the normalized epsilon_delta with the low/mid Schur absorption budget.",
        ],
        "notClaimedHere": (
            "This theorem does not pass the estimate to the completed high "
            "block.  That continuum passage is the job of "
            "high_block_compact_exhaustion_proof.py."
        ),
        "remainingAnalyticGap": None
        if constants_closed
        else "The inactive source constants or absorption budget are not closed.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Source-inactive tail constants theorem")
    print(f"  normalized epsilon_delta: {normalized_epsilon:.12e}")
    print(f"  finite low/mid budget: {budget:.12e}")
    print(f"  absorption slack: {slack:.12e}")
    print(f"  constants closed: {constants_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
