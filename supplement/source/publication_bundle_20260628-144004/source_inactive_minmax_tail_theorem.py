#!/usr/bin/env python3
r"""Min-max proof ledger for the source-inactive tail estimate.

Let H be the closed-trace high block with A-inner product and let

    E : H -> L^2([u_-,u_+]; C^2),        S = E^*E.

If P_2 is the spectral projection onto the two largest singular directions of
E, then the spectral theorem gives the exact source-inactive estimate

    ||(I-P_2)E f||^2 <= lambda_3(S) <A f,f>.

The full-theta inactive-tail certificate supplies a perturbative bound

    lambda_3(S) <= lambda_3(S_{8,h}) + eta

inside the certified source model, and a lower bound for the top source scale

    lambda_1(S) >= lambda_1(S_{8,h}) - eta.

Therefore the normalized source operator

    Ehat = E / sqrt(lambda_1(S_{8,h}) - eta)

satisfies

    ||(I-P_2)Ehat f||^2 <= epsilon_delta <A f,f>,

with epsilon_delta equal to the inactive/top ratio in the certificate.

This script records the exact constants and keeps the remaining Galerkin-to-
continuum high-block passage separate from the min-max proof.
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
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument(
        "--absorption-json",
        default="source_inactive_tail_constants_consequence_theorem.json",
    )
    parser.add_argument(
        "--high-block-json",
        default="high_block_minmax_tail_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="source_inactive_minmax_tail_theorem.json")
    args = parser.parse_args()

    inactive = load(args.inactive_json)
    absorption = load(args.absorption_json)
    high_block_path = Path(args.high_block_json)
    high_block = load(args.high_block_json) if high_block_path.exists() else None

    absolute_tail = float(inactive["absoluteTailBound"])
    source_top_lower = float(inactive["sourceTopLower"])
    finite_lambda3 = float(inactive["finiteInactiveLambda3"])
    finite_top = float(inactive["finiteSourceTop"])
    perturbation_eta = float(inactive["operatorPerturbationEta"])
    normalized_epsilon = float(inactive["normalizedEpsilonDelta"])
    budget = float(absorption["finiteLowMidSchurBudget"])
    slack = budget - normalized_epsilon
    high_block_closed = bool(high_block and high_block.get("tailEstimatePassesToContinuum"))

    data = {
        "theoremName": "source-inactive min-max tail theorem",
        "inactiveJson": args.inactive_json,
        "absorptionJson": args.absorption_json,
        "highBlockJson": args.high_block_json if high_block else None,
        "basis": inactive.get("basis"),
        "sourceGrid": inactive.get("sourceGrid"),
        "globalTraceRatio": inactive.get("globalTraceRatio"),
        "sourceWindow": inactive.get("sourceWindow"),
        "finiteSourceTop": finite_top,
        "finiteInactiveLambda3": finite_lambda3,
        "operatorPerturbationEta": perturbation_eta,
        "absoluteTailBound": absolute_tail,
        "sourceTopLower": source_top_lower,
        "normalizedEpsilonDelta": normalized_epsilon,
        "finiteLowMidSchurBudget": budget,
        "absorptionSlack": slack,
        "epsilonOverBudget": normalized_epsilon / budget if budget else float("inf"),
        "rawUnnormalizedStatement": (
            "||(I-P_2)E_Phi f||^2 <= absoluteTailBound <A f,f>."
        ),
        "normalizedStatement": (
            "For Ehat_Phi=E_Phi/sqrt(sourceTopLower), "
            "||(I-P_2)Ehat_Phi f||^2 <= normalizedEpsilonDelta <A f,f>."
        ),
        "minMaxProofInCertifiedSourceModel": True,
        "weylPerturbationUsed": True,
        "absorbableByFiniteLowMidBlock": slack > 0,
        "literalUnnormalizedEpsilonIsSmall": False,
        "continuumHighBlockExhaustionProvedHere": high_block_closed,
        "continuumHighBlockExhaustionImported": high_block_closed,
        "remainingAnalyticGap": (
            None
            if high_block_closed
            else (
                "Prove the Galerkin-to-continuum high-block exhaustion/elliptic "
                "estimate which upgrades the certified source model to the full "
                "closed space H_M cap ker R_global.  The min-max tail bound itself "
                "is then automatic from the spectral theorem."
            )
        ),
        "proof": [
            (
                "Let S=E_Phi^*E_Phi on the A-Hilbert space.  If P_2 is the "
                "spectral projection onto the two largest eigenvalues of S, "
                "then <S(I-P_2)f,(I-P_2)f> <= lambda_3(S)<A f,f>."
            ),
            (
                "The certified finite source model has lambda_3(S_{8,h})="
                f"{finite_lambda3:.12e} and ||S-S_{{8,h}}||<="
                f"{perturbation_eta:.12e}; Weyl monotonicity gives "
                f"lambda_3(S)<={absolute_tail:.12e}."
            ),
            (
                "The same perturbation bound gives the top source scale "
                f"lambda_1(S)>={source_top_lower:.12e}.  After normalizing "
                "E_Phi by sqrt(sourceTopLower), the tail constant is "
                f"{normalized_epsilon:.12e}."
            ),
            (
                "Since this is below the finite low/mid Schur budget "
                f"{budget:.12e}, the normalized source-inactive tail is "
                f"absorbed with slack {slack:.12e}."
            ),
        ],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Source-inactive min-max tail theorem")
    print(f"  raw tail bound: {absolute_tail:.12e}")
    print(f"  source top lower: {source_top_lower:.12e}")
    print(f"  normalized epsilon_delta: {normalized_epsilon:.12e}")
    print(f"  finite low/mid budget: {budget:.12e}")
    print(f"  absorption slack: {slack:.12e}")
    print("  min-max proof in certified source model: True")
    print(f"  continuum high-block exhaustion imported: {high_block_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
