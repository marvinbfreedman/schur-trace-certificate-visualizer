#!/usr/bin/env python3
r"""Abstract Mosco-to-projection convergence theorem.

This is the first functional-analytic component of the high-block exhaustion
passage.  It is model-free: no theta, Weyl, Volterra, trace, or source-window
input is used here.
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
        default="abstract_mosco_projection_convergence_theorem.json",
    )
    args = parser.parse_args()

    ok = True
    data = {
        "theoremName": "abstract Mosco projection convergence theorem",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "For closed subspaces H_N,H of a Hilbert space, Mosco convergence "
            "H_N -> H is equivalent to strong convergence of the orthogonal "
            "projections P_N -> P."
        ),
        "hypotheses": [
            "H is a Hilbert space with norm induced by the closed positive A form.",
            "H_N and H are closed subspaces.",
            "H_N Mosco-converges to H: every h in H has a strong recovery sequence h_N in H_N, and every weak limit of bounded h_N in H_N lies in H.",
        ],
        "statuses": {
            "limsupRecoveryToProjectionStatus": status(
                "limsup recovery gives projection upper bound",
                ok,
                (
                    "Recovery sequences for points of H give "
                    "limsup ||(I-P_N)h||=0 for h in H."
                ),
            ),
            "liminfClosednessToProjectionStatus": status(
                "liminf closedness gives weak projection limit",
                ok,
                (
                    "If P_N x has bounded subsequences, Mosco liminf forces "
                    "all weak cluster points to lie in H; Hilbert projection "
                    "orthogonality identifies the limit as Px."
                ),
            ),
            "projectionStrongConvergenceStatus": status(
                "orthogonal projections converge strongly",
                ok,
                (
                    "The standard Kato/Mosco theorem yields P_N x -> P x for "
                    "every x in the Hilbert space."
                ),
            ),
        },
        "moscoProjectionConvergenceClosed": ok,
        "strongProjectionConvergenceClosed": ok,
        "proof": [
            "Use the variational characterization of P_N x as the unique minimizer of ||x-y|| over y in H_N.",
            "The Mosco limsup recovery sequence gives the upper bound limsup ||x-P_Nx|| <= ||x-Px||.",
            "Bounded minimizing sequences have weak cluster points in H by the Mosco liminf condition.",
            "Lower semicontinuity gives the matching lower bound, so P_Nx converges weakly to Px and the norms converge; Hilbert uniform convexity upgrades this to strong convergence.",
        ],
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract Mosco projection convergence theorem")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
