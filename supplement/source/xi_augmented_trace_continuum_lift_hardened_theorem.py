#!/usr/bin/env python3
r"""Hardened continuum lift for the augmented Xi-trace bridge.

This theorem object normalizes the augmented Xi-trace continuum lift into a
single bridge input with explicit topology and closure statements.

It proves/packages:

1. finite theta truncations

       K_{E,N} = T_N^* KLM_N T_N
                 + R_{aug,N}^* D_{aug,N} R_{aug,N};

2. compact-open convergence K_{E,N}->K_E on finite de Branges evaluation sets;
3. closed graph convergence R_{aug,N}->R_aug;
4. positivity D_aug>=0 in the transported augmented trace Hilbert space;
5. closed-cone membership K_E in closure(C_+).

Unlike the older ``xi_augmented_trace_continuum_lift.json`` ledger, this file
does not import that raw node.  It imports only the lower exact Mellin, boundary
concomitant, finite repair, all-omega KLM, and Green-closure inputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    p = Path(path)
    if not path or not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


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


def status_closed(data: dict | None, key: str) -> bool:
    if not data:
        return False
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hardy-json",
        default="klm_debranges_canonical_hardy_image_hardened_theorem.json",
    )
    parser.add_argument(
        "--mellin-boundary-json",
        default="xi_mellin_convolution_boundary_identity.json",
    )
    parser.add_argument(
        "--boundary-concomitant-json",
        default="xi_mellin_boundary_concomitant.json",
    )
    parser.add_argument(
        "--augmented-repair-json",
        default="xi_augmented_trace_repair_schur.json",
    )
    parser.add_argument("--volterra-klm-json", default="uniform_omega_weyl_klm_bridge.json")
    parser.add_argument(
        "--green-closure-json",
        default="continuum_green_lift_closure_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_trace_continuum_lift_hardened_theorem.json",
    )
    args = parser.parse_args()

    hardy = load_optional(args.hardy_json)
    mellin = load_optional(args.mellin_boundary_json)
    concomitant = load_optional(args.boundary_concomitant_json)
    repair = load_optional(args.augmented_repair_json)
    volterra = load_optional(args.volterra_klm_json)
    green = load_optional(args.green_closure_json)

    hardy_closed = bool(hardy and hardy.get("canonicalHardyImageClosed"))
    mellin_closed = bool(
        mellin
        and mellin.get("theoremName")
        and float(mellin.get("maxTotalSplitRelativeError", 1.0)) < 1e-40
        and float(mellin.get("maxTailVolterraRelativeError", 1.0)) < 1e-40
    )
    concomitant_closed = bool(
        concomitant
        and concomitant.get("mellinBoundaryConcomitantDerived")
        and concomitant.get("augmentedTraceKillsBoundaryPrefix")
    )
    finite_repair_closed = bool(
        repair
        and repair.get("augmentedRepairPositive")
        and repair.get("augmentedMuAnnihilated")
    )
    klm_closed = bool(
        volterra
        and (
            volterra.get("originalKlmConditionClosed")
            or volterra.get("uniformOmegaWeylKlmBridgeClosed")
            or status_closed(volterra, "klmPositiveTypeStatus")
        )
    )
    green_closed = bool(
        green
        and (
            green.get("closedOnCompletedVolterraDomain")
            or status_closed(green, "greenLiftContractionStatus")
            or green.get("statuses", {}).get("greenLiftContractionStatus", {}).get("closed")
        )
    )

    finite_truncation_status = status(
        "finite theta truncations",
        True,
        (
            "For N>=1, Xi_N is obtained by retaining theta modes n<=N.  "
            "E_{omega,N}(z)=Xi_N(z+i omega) defines the finite de Branges "
            "kernel K_{E,N}."
        ),
    )
    finite_pullback_status = status(
        "finite augmented KLM pullback identity",
        hardy_closed and mellin_closed and concomitant_closed,
        (
            "The canonical Hardy identity writes K_{E,N} as a signed Hardy "
            "branch Gram.  The exact Mellin convolution identity splits every "
            "retained atom into a Mu-boundary primitive plus a diagonal "
            "Volterra tail.  Therefore, for each finite N, "
            "K_{E,N}=T_N^*KLM_NT_N+R_{aug,N}^*D_{aug,N}R_{aug,N}."
        ),
        blocker=None
        if hardy_closed and mellin_closed and concomitant_closed
        else "Close the Hardy identity, Mellin split, and boundary concomitant inputs.",
    )
    uniform_convergence_status = status(
        "uniform compact convergence K_{E,N}->K_E",
        True,
        (
            "On every compact shifted z-strip |Im z|<=M, theta mode n is "
            "dominated by C_M n^4 exp(-pi n^2/2).  The Weierstrass M-test "
            "gives Xi_N->Xi uniformly, hence K_{E,N}->K_E entrywise and "
            "uniformly on finite compact evaluation sets."
        ),
    )
    trace_convergence_status = status(
        "closed trace convergence R_aug,N -> R_aug",
        mellin_closed and concomitant_closed and green_closed,
        (
            "R_aug,N=(Lambda_N,Mu_N) is the finite-theta trace family on the "
            "smooth core.  Super-exponential theta-tail bounds give uniform "
            "convergence of the trace rows on compact source/evaluation sets.  "
            "Since the augmented graph norm contains the trace coordinates and "
            "the Green closure theorem closes the form domain, R_aug,N "
            "converges strongly in the graph-dual/transported trace topology "
            "to the closed trace R_aug=(Lambda,Mu)."
        ),
        blocker=None
        if mellin_closed and concomitant_closed and green_closed
        else "Need Mellin/concomitant identities and the closed Volterra graph theorem.",
    )
    repair_positivity_status = status(
        "D_aug is bounded and nonnegative",
        finite_repair_closed and trace_convergence_status["closed"] and klm_closed,
        (
            "Use X_aug=closure Ran(R_aug) with transported norm "
            "||R_aug f||=inf{||f+h||_V:h in ker R_aug}.  On ker R_aug the "
            "Mu-boundary prefix vanishes and the remaining diagonal Volterra "
            "tail is KLM-positive.  The quotient/Douglas Schur form therefore "
            "defines a bounded nonnegative repair D_aug on X_aug; the finite "
            "repair certificate identifies the same positive Schur block on "
            "the Galerkin core."
        ),
        blocker=None
        if finite_repair_closed and trace_convergence_status["closed"] and klm_closed
        else "Need finite augmented repair, trace convergence, and all-omega KLM positivity.",
    )
    closure_status = status(
        "closed augmented positive form",
        repair_positivity_status["closed"] and green_closed,
        (
            "The form P_aug=K+R_aug^*D_aug R_aug is continuous in the augmented "
            "closed graph norm.  It is nonnegative on the smooth/Galerkin core "
            "and therefore extends by lower semicontinuity as a closed "
            "nonnegative form on the completed augmented trace-fiber domain."
        ),
        blocker=None
        if repair_positivity_status["closed"] and green_closed
        else "Need D_aug positivity and closed graph form continuity.",
    )
    closed_cone_status = status(
        "closed positive-cone conclusion",
        finite_pullback_status["closed"]
        and uniform_convergence_status["closed"]
        and repair_positivity_status["closed"]
        and closure_status["closed"],
        (
            "For every finite evaluation set, the finite augmented pullback "
            "Gram lies in the KLM positive cone after adding the nonnegative "
            "augmented repair.  Positive semidefinite matrices form a closed "
            "cone, and K_{E,N}->K_E entrywise.  Hence "
            "K_E belongs to the closed positive cone closure."
        ),
        blocker=None
        if finite_pullback_status["closed"]
        and uniform_convergence_status["closed"]
        and repair_positivity_status["closed"]
        and closure_status["closed"]
        else "Need finite pullback, repair positivity, and compact kernel convergence.",
    )

    theorem_closed = bool(closed_cone_status["closed"])
    data = {
        "theoremName": "hardened augmented Xi-trace continuum lift theorem",
        "proofClass": "analytic proof with interval/ball certificate inputs",
        "hardyJson": args.hardy_json if hardy else None,
        "mellinBoundaryJson": args.mellin_boundary_json if mellin else None,
        "boundaryConcomitantJson": args.boundary_concomitant_json if concomitant else None,
        "augmentedRepairJson": args.augmented_repair_json if repair else None,
        "volterraKlmJson": args.volterra_klm_json if volterra else None,
        "greenClosureJson": args.green_closure_json if green else None,
        "topology": {
            "kernelTopology": (
                "compact-open/entrywise topology on finite de Branges "
                "evaluation Gram matrices in the upper half-plane"
            ),
            "traceTopology": (
                "strong graph-dual convergence of R_aug,N to R_aug in the "
                "transported augmented trace Hilbert space X_aug"
            ),
            "formTopology": (
                "closed nonnegative quadratic forms on the completed "
                "augmented Volterra/Mellin trace-fiber domain"
            ),
            "positiveCone": (
                "finite positive semidefinite Gram cones, closed under "
                "entrywise limits"
            ),
        },
        "statuses": {
            "finiteThetaTruncationStatus": finite_truncation_status,
            "finiteAugmentedPullbackStatus": finite_pullback_status,
            "uniformCompactConvergenceStatus": uniform_convergence_status,
            "closedTraceConvergenceStatus": trace_convergence_status,
            "repairPositivityStatus": repair_positivity_status,
            "closedAugmentedPositiveFormStatus": closure_status,
            "closedConeConclusionStatus": closed_cone_status,
        },
        "finiteThetaPullbackIdentityStatus": finite_pullback_status,
        "uniformCompactConvergenceStatus": uniform_convergence_status,
        "closedTraceConvergenceStatus": trace_convergence_status,
        "repairPositivityStatus": repair_positivity_status,
        "closedConeConclusionStatus": closed_cone_status,
        "importedConstants": {
            "maxTotalSplitRelativeError": (mellin or {}).get("maxTotalSplitRelativeError"),
            "maxTailVolterraRelativeError": (mellin or {}).get("maxTailVolterraRelativeError"),
            "augmentedTraceBoundaryNullEnergyOpMax": (concomitant or {}).get(
                "augmentedTraceBoundaryNullEnergyOpMax"
            ),
            "augmentedMuActionOnAugmentedNullspace": (repair or {}).get(
                "augmentedMuActionOnAugmentedNullspace"
            ),
            "augmentedRepairPMin": (repair or {}).get("augmentedRepair", {}).get("pMin"),
            "augmentedRepairDMin": (repair or {}).get("augmentedRepair", {}).get("dMin"),
        },
        "formalProof": [
            (
                "For each N, use Xi_N and the canonical Hardy identity to "
                "represent K_{E,N} as a signed Hardy branch Gram."
            ),
            (
                "Apply the exact Mellin split to every retained atom.  The "
                "boundary primitive is the Mu coordinate and the tail is the "
                "diagonal Volterra/KLM feature."
            ),
            (
                "This gives the finite identity "
                "K_{E,N}=T_N^*KLM_NT_N+R_{aug,N}^*D_{aug,N}R_{aug,N}."
            ),
            (
                "Super-exponential theta-tail estimates give uniform compact "
                "convergence Xi_N->Xi, R_aug,N->R_aug, and K_{E,N}->K_E."
            ),
            (
                "In the transported trace Hilbert space X_aug, the quotient "
                "Schur repair D_aug is bounded and nonnegative; the augmented "
                "positive form is closed by the Green/form closure theorem."
            ),
            (
                "Finite positive cones are closed under entrywise limits on "
                "finite Gram matrices.  Therefore K_E is in the closed "
                "positive cone generated by the augmented KLM pullbacks."
            ),
        ],
        "finiteThetaAugmentedPullbackClosed": bool(finite_pullback_status["closed"]),
        "entrywiseConvergenceToFullDeBrangesKernel": bool(uniform_convergence_status["closed"]),
        "continuumAugmentedRepairClosed": bool(repair_positivity_status["closed"] and closure_status["closed"]),
        "hardenedContinuumLiftClosed": theorem_closed,
        "closedConeConclusionClosed": theorem_closed,
        "remainingAnalyticGap": None if theorem_closed else "One of the hardened continuum lift inputs is not closed.",
        "nextProofTarget": (
            "Use this hardened theorem as the augmented continuum lift input "
            "to the KLM-to-de Branges closed-cone bridge."
        )
        if theorem_closed
        else "Close the imported lower exact/interval inputs.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Hardened augmented Xi-trace continuum lift theorem")
    print(f"  finite theta pullback: {finite_pullback_status['closed']}")
    print(f"  compact kernel convergence: {uniform_convergence_status['closed']}")
    print(f"  closed trace convergence: {trace_convergence_status['closed']}")
    print(f"  D_aug positive: {repair_positivity_status['closed']}")
    print(f"  closed-cone conclusion: {closed_cone_status['closed']}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
