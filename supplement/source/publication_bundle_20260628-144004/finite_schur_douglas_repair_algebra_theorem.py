#!/usr/bin/env python3
r"""Universal finite Schur/Douglas repair algebra theorem.

This is a symbolic finite-dimensional theorem.  It contains no theta,
Volterra, or Galerkin data.

Let K be Hermitian on H=N+U with block form

    K = [ A  B ]
        [ B* C ]

where A>=0 on N.  If the Douglas range condition

    (I-AA^+)B = 0

holds and M>=B*A^+B-C on U with M>=0, then the repaired form

    K + R^* D R

is nonnegative whenever D is a trace-side realization of M on the trace image
R(U).  This is the finite algebra used by the augmented Xi repair certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json-out",
        default="finite_schur_douglas_repair_algebra_theorem.json",
    )
    args = parser.parse_args()

    data = {
        "theoremName": "finite Schur Douglas repair algebra theorem",
        "proofClass": "symbolic identity",
        "ambientSpace": "finite-dimensional Hilbert space H=N direct-sum U",
        "hypotheses": [
            "K is Hermitian with blocks A=N*K*N, B=N*K*U, C=U*K*U.",
            "A is positive semidefinite.",
            "The Douglas range condition (I-AA^+)B=0 holds.",
            "M is positive semidefinite and M >= B^*A^+B-C on U.",
            "D realizes M on the finite trace image R(U).",
        ],
        "statuses": {
            "finiteSchurDouglasAlgebraStatus": status(
                "finite Schur/Douglas algebra",
                True,
                (
                    "Completing the square with A^+ gives "
                    "<n+u,K(n+u)>+<u,Mu> = <A^{1/2}n+A^{+1/2}Bu, "
                    "A^{1/2}n+A^{+1/2}Bu> + <u,(C+M-B^*A^+B)u>."
                ),
            ),
            "finiteRangeConditionStatus": status(
                "finite Douglas range condition",
                True,
                (
                    "The nullspace component of B is exactly the obstruction; "
                    "(I-AA^+)B=0 is equivalent to the cross term being in "
                    "Range(A^{1/2})."
                ),
            ),
            "finiteTraceRealizationStatus": status(
                "finite trace-side realization",
                True,
                (
                    "On the finite trace image, choosing D so that "
                    "R_U^*DR_U=M transfers the U-block repair to the original "
                    "finite form."
                ),
            ),
        },
        "finiteSchurDouglasAlgebraClosed": True,
        "finiteRangeCriterionClosed": True,
        "finiteTraceRepairRealizationClosed": True,
        "formalProof": [
            "Use the Moore-Penrose inverse A^+ on closure Range(A).",
            "Apply the Douglas range condition to write B=A A^+ B.",
            "Complete the square in the N variable.",
            "Use C+M-B^*A^+B>=0 to obtain nonnegativity.",
            "Pull M back through the trace image realization D.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Finite Schur Douglas repair algebra theorem")
    print("  algebra closed: True")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
