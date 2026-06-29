#!/usr/bin/env python3
r"""Repair-free Schur equivalence ledger.

Let the transported Volterra form be blocked relative to

    V = N \oplus U,      N = ker R_global,

as

    q(n+u) = a(n) + 2 Re b(n,u) + c(u).

Assume the closed quotient hypotheses already proved:

    a >= 0 on N,
    b(n,u)=<A^{1/2}n, Gamma u>.

Then

    q(n+u)
      = ||A^{1/2}n + Gamma u||^2
        + <(C-Gamma^*Gamma)u,u>.

Therefore the repair-free theorem

    D_q = (Gamma^*Gamma-C)_+ = 0

is equivalent to

    Gamma^*Gamma <= C,

which is equivalent to q>=0 on the whole completed form domain.  This is the
direct full-form positivity theorem; no boundary repair is left after the
primitive transport audit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quotient-json",
        default="weyl_volterra_quotient_schur_theorem.json",
    )
    parser.add_argument(
        "--repair-summary-json",
        default="dq_vanishing_repair_route_summary.json",
    )
    parser.add_argument(
        "--primitive-endpoint-json",
        default="primitive_endpoint_compatibility_theorem.json",
    )
    parser.add_argument("--json-out", default="repair_free_schur_equivalence.json")
    args = parser.parse_args()

    quotient = load(args.quotient_json)
    repair = load(args.repair_summary_json)
    primitive_endpoint = load(args.primitive_endpoint_json)

    positivity_on_n = bool(quotient.get("positivityOnKerRStatus", {}).get("closed"))
    douglas = bool(quotient.get("douglasCrossFormStatus", {}).get("closed"))
    primitive_zero = bool(
        repair.get("statuses", {})
        .get("primitiveBoundaryZeroStatus", {})
        .get("closed")
    )
    finite_support = bool(
        repair.get("statuses", {})
        .get("finiteSchurDefectStatus", {})
        .get("closed")
    )
    endpoint_dq_zero = bool(
        primitive_endpoint.get("statuses", {})
        .get("dqVanishesOnXRStatus", {})
        .get("closed")
    )

    quotient_hypotheses = status(
        "quotient Schur hypotheses for block factorization",
        positivity_on_n and douglas,
        (
            "The closed quotient theorem gives a>=0 on N=ker R and the "
            "Douglas factorization b(n,u)=<A^{1/2}n,Gamma u>."
        ),
    )
    algebraic_equivalence = status(
        "D_q=0 iff full transported form is positive",
        quotient_hypotheses["closed"],
        (
            "Completing the square gives "
            "q(n+u)=||A^{1/2}n+Gamma u||^2+<"
            "(C-Gamma^*Gamma)u,u>.  Hence q>=0 on N+U iff "
            "Gamma^*Gamma<=C iff D_q=(Gamma^*Gamma-C)_+=0."
        ),
        blocker=None
        if quotient_hypotheses["closed"]
        else "Close positivity on ker R and the Douglas factorization.",
    )
    primitive_lift = status(
        "primitive lift needs no boundary repair",
        primitive_zero,
        (
            "The primitive transport audit gives D_bdy=0, so original "
            "positivity follows from the repair-free full-form positivity "
            "theorem without an added endpoint term."
        ),
    )
    finite_evidence = status(
        "finite repair-free evidence",
        finite_support,
        (
            "All imported finite scans have lambda_max(Gamma^*Gamma-C)<0."
        ),
    )
    continuum = status(
        "continuum repair-free theorem",
        endpoint_dq_zero,
        (
            (
                "The primitive endpoint compatibility theorem identifies the "
                "Schur trace form as P-M and uses the closed Green-lift "
                "contraction M<=P.  Hence Gamma^*Gamma<=C and D_q=0 on X_R."
            )
            if endpoint_dq_zero
            else (
                "The algebraic equivalence is closed and finite evidence is "
                "negative, but the continuum inequality Gamma^*Gamma<=C on "
                "X_R has not been proved."
            )
        ),
        blocker=None
        if endpoint_dq_zero
        else (
            "Prove C-Gamma^*Gamma>=0 on X_R, preferably by an explicit "
            "Volterra/Gram factorization of the Schur complement."
        ),
    )

    data = {
        "theoremName": "repair-free Schur equivalence",
        "quotientJson": args.quotient_json,
        "repairSummaryJson": args.repair_summary_json,
        "primitiveEndpointJson": args.primitive_endpoint_json,
        "statuses": {
            "quotientHypothesesStatus": quotient_hypotheses,
            "algebraicEquivalenceStatus": algebraic_equivalence,
            "primitiveLiftNoBoundaryRepairStatus": primitive_lift,
            "finiteRepairFreeEvidenceStatus": finite_evidence,
            "continuumRepairFreeTheoremStatus": continuum,
        },
        "identity": (
            "q(n+u)=||A^{1/2}n+Gamma u||^2+<"
            "(C-Gamma^*Gamma)u,u>"
        ),
        "equivalentTargets": [
            "D_q=0 on X_R",
            "Gamma^*Gamma<=C on X_R",
            "C-Gamma^*Gamma>=0 on X_R",
            "Q_Phi>=0 on the full completed primitive/form domain",
        ],
        "finiteEvidence": {
            "scanRows": repair.get("scanRows", []),
            "worstFiniteRow": repair.get("worstFiniteRow"),
            "allFiniteDqZero": repair.get("allFiniteDqZero"),
        },
        "continuumTheoremClosed": continuum["closed"],
        "nextProofTarget": continuum["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Repair-free Schur equivalence")
    print(f"  quotient hypotheses: {quotient_hypotheses['closed']}")
    print(f"  algebraic equivalence: {algebraic_equivalence['closed']}")
    print(f"  primitive D_bdy=0: {primitive_lift['closed']}")
    print(f"  finite evidence: {finite_evidence['closed']}")
    print(f"  continuum theorem: {continuum['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
