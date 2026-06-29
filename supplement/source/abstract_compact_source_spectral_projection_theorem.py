#!/usr/bin/env python3
r"""Abstract Riesz spectral projection stability theorem for compact sources."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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
        "--compression-json",
        default="abstract_compact_compression_norm_consequence_theorem.json",
    )
    parser.add_argument(
        "--riesz-json",
        default="abstract_riesz_projection_continuity_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_compact_source_spectral_projection_theorem.json",
    )
    args = parser.parse_args()

    compression = load(args.compression_json)
    riesz = load(args.riesz_json)
    compression_ok = bool(compression.get("compactCompressionNormConvergenceClosed"))
    riesz_ok = bool(riesz.get("rieszProjectionContinuityClosed"))
    ok = compression_ok and riesz_ok

    data = {
        "theoremName": "abstract compact source spectral projection theorem",
        "proofClass": "abstract functional analytic lemma",
        "statement": (
            "For compact self-adjoint positive operators S_N -> S in norm, "
            "isolated spectral clusters have convergent Riesz projections and "
            "stable finite-rank active/inactive splittings."
        ),
        "hypotheses": [
            "S_N and S are compact self-adjoint positive operators.",
            "||S_N-S|| -> 0.",
            "The selected source cluster is separated from the complement by a positive contour gap.",
        ],
        "statuses": {
            "normConvergenceInputStatus": status(
                "compact-source norm convergence input",
                compression_ok,
                "The compact compression theorem gives norm convergence after subspace compression.",
            ),
            "rieszProjectionStabilityStatus": status(
                "Riesz spectral projection stability",
                riesz_ok,
                (
                    "The abstract Riesz projection continuity theorem gives "
                    "norm convergence of the spectral projections once a "
                    "contour gap is supplied."
                ),
            ),
            "activeInactiveSplittingStabilityStatus": status(
                "active/inactive splitting stability",
                ok,
                (
                    "The finite active rank and inactive orthogonal complement "
                    "are stable under norm-small compact self-adjoint perturbations."
                ),
            ),
        },
        "compactSourceRieszProjectionConvergenceClosed": ok,
        "sourceActiveInactiveSplittingStabilityClosed": ok,
        "proof": [
            "Choose a contour enclosing the isolated source cluster and no other spectral point of S.",
            "For N large, ||S_N-S|| is smaller than the contour distance, so the resolvents exist on the same contour.",
            "Integrate the resolvent identity to obtain Riesz projection norm convergence.",
            "The rank is constant for large N and the active/inactive projections converge in norm.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Compact-source norm convergence or Riesz projection input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract compact source spectral projection theorem")
    print(f"  compression input: {compression_ok}")
    print(f"  Riesz input: {riesz_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
