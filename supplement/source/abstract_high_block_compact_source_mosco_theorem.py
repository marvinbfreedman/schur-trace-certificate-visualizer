#!/usr/bin/env python3
r"""Abstract high-block Mosco/compact-source theorem.

This is a model-free Hilbert-space lemma.  It does not assert any theta,
Volterra, source-window, or finite certificate estimate.

If finite closed subspaces H_N Mosco-converge to H in an A-Hilbert space and
S is a compact positive source operator, then the A-orthogonal projections
Pi_N converge strongly to Pi and

    ||Pi_N S Pi_N - Pi S Pi||_{A->A} -> 0.

Consequently spectral min-max tail inequalities for the compact positive
operators pass from the finite H_N model to the continuum H model.
"""

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
        "--projection-json",
        default="abstract_mosco_projection_convergence_theorem.json",
    )
    parser.add_argument(
        "--compression-json",
        default="abstract_compact_compression_norm_convergence_theorem.json",
    )
    parser.add_argument(
        "--spectral-json",
        default="abstract_compact_source_spectral_projection_theorem.json",
    )
    parser.add_argument(
        "--minmax-json",
        default="abstract_minmax_tail_passage_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="abstract_high_block_compact_source_mosco_theorem.json",
    )
    args = parser.parse_args()

    projection = load(args.projection_json)
    compression = load(args.compression_json)
    spectral = load(args.spectral_json)
    minmax = load(args.minmax_json)
    projection_ok = bool(projection.get("moscoProjectionConvergenceClosed"))
    compression_ok = bool(compression.get("compactCompressionNormConvergenceClosed"))
    spectral_ok = bool(spectral.get("compactSourceRieszProjectionConvergenceClosed"))
    minmax_ok = bool(minmax.get("abstractMinmaxTailPassageClosed"))
    ok = all([projection_ok, compression_ok, spectral_ok, minmax_ok])
    data = {
        "theoremName": "abstract high-block compact-source Mosco theorem",
        "proofClass": "symbolic identity",
        "projectionJson": args.projection_json,
        "compressionJson": args.compression_json,
        "spectralJson": args.spectral_json,
        "minmaxJson": args.minmax_json,
        "statement": (
            "Mosco convergence of closed A-Hilbert subspaces plus compactness "
            "of the source operator implies source-operator norm convergence, "
            "Riesz spectral projection convergence, and min-max tail passage."
        ),
        "statuses": {
            "moscoProjectionConvergencePrincipleStatus": status(
                "Mosco projection convergence principle",
                projection_ok,
                (
                    "Imported from the abstract Mosco projection convergence "
                    "theorem."
                ),
            ),
            "compactSourceNormConvergencePrincipleStatus": status(
                "compact source norm convergence principle",
                compression_ok,
                (
                    "Imported from the abstract compact-compression norm "
                    "convergence theorem."
                ),
            ),
            "spectralMinmaxPassagePrincipleStatus": status(
                "spectral min-max passage principle",
                spectral_ok and minmax_ok,
                (
                    "Imported from the abstract compact-source spectral "
                    "projection theorem and abstract min-max tail passage "
                    "theorem."
                ),
            ),
        },
        "abstractHighBlockCompactSourceMoscoClosed": ok,
        "moscoProjectionConvergenceClosed": projection_ok,
        "compactSourceNormConvergenceClosed": compression_ok,
        "compactSourceRieszProjectionConvergenceClosed": spectral_ok,
        "spectralMinmaxPassageClosed": minmax_ok,
        "noModelSpecificInputClaimed": True,
        "proof": [
            "This file is now a legacy wrapper over the split abstract functional-analysis lemmas.",
            "The Mosco-to-projection, compact-compression, Riesz projection, and min-max steps are audited separately.",
        ],
        "remainingAnalyticGap": None if ok else "One split abstract functional-analysis lemma is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Abstract high-block compact-source Mosco theorem")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
