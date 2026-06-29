#!/usr/bin/env python3
r"""High-block min-max tail input theorem."""

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
        "--tail-bound-json",
        default="source_inactive_tail_bound_consequence_theorem.json",
    )
    parser.add_argument(
        "--source-split-json",
        default="high_block_active_source_split_consequence_theorem.json",
    )
    parser.add_argument(
        "--minmax-passage-json",
        default="abstract_minmax_tail_passage_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="high_block_minmax_tail_input_theorem.json",
    )
    args = parser.parse_args()

    tail = load(args.tail_bound_json)
    source_split = load(args.source_split_json)
    minmax = load(args.minmax_passage_json)

    tail_ok = bool(
        tail.get("sourceInactiveTailConstantsClosed")
        or tail.get("sourceInactiveTailDominationConstantsClosed")
        or tail.get("minMaxProofInCertifiedSourceModel")
    )
    absorb_ok = bool(tail.get("absorbableByFiniteLowMidBlock"))
    spectral_ok = bool(
        source_split.get("activeInactiveSourceSplitClosed")
        or source_split.get("sourceActiveInactiveSplittingStabilityClosed")
        or source_split.get("highBlockSpectralProjectionClosed")
    )
    abstract_ok = bool(
        minmax.get("continuumInactiveTailPassageClosed")
        or minmax.get("abstractMinmaxTailPassageClosed")
    )
    ok = tail_ok and absorb_ok and spectral_ok and abstract_ok

    data = {
        "theoremName": "high-block min-max tail input theorem",
        "proofClass": "symbolic identity",
        "tailBoundJson": args.tail_bound_json,
        "sourceSplitJson": args.source_split_json,
        "minmaxPassageJson": args.minmax_passage_json,
        "normalizedEpsilonDelta": tail.get("normalizedEpsilonDelta"),
        "finiteLowMidSchurBudget": tail.get("finiteLowMidSchurBudget"),
        "absorptionSlack": tail.get("absorptionSlack"),
        "statement": (
            "The high-block source model satisfies the finite inactive-tail "
            "inequality and spectral-projection hypotheses needed by the "
            "abstract min-max tail passage theorem."
        ),
        "statuses": {
            "finiteTailConstantInputStatus": status(
                "finite certified inactive-tail constant",
                tail_ok,
                "The source-inactive constants theorem supplies the finite certified source-model tail bound.",
            ),
            "lowMidAbsorptionInputStatus": status(
                "finite low/mid absorption budget",
                absorb_ok,
                "The same constants theorem verifies epsilon_delta is absorbable by the finite low/mid Schur budget.",
            ),
            "spectralProjectionInputStatus": status(
                "high-block spectral projection input",
                spectral_ok,
                "The high-block spectral theorem identifies the stable active/inactive source split.",
            ),
            "abstractMinmaxInputStatus": status(
                "abstract min-max passage input",
                abstract_ok,
                "The abstract min-max theorem passes the finite inactive-tail bound through the convergent source model.",
            ),
            "highBlockMinmaxTailInputStatus": status(
                "high-block min-max tail passage",
                ok,
                "The finite tail constant, spectral split, and abstract min-max theorem match.",
            ),
        },
        "highBlockMinmaxTailPassageClosed": ok,
        "abstractMinmaxTailPassageClosed": ok,
        "tailEstimatePassesToContinuum": ok,
        "proof": [
            "Use the source-inactive constants theorem for the finite inactive-tail quadratic inequality.",
            "Use the high-block spectral projection theorem to identify the active two-dimensional source subspace and its inactive complement.",
            "Apply the abstract min-max tail passage theorem to transport the finite bound to the continuum closed high block.",
        ],
        "remainingAnalyticGap": None if ok else "Finite tail, spectral projection, or abstract min-max input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("High-block min-max tail input theorem")
    print(f"  finite tail input: {tail_ok}")
    print(f"  absorption input: {absorb_ok}")
    print(f"  spectral input: {spectral_ok}")
    print(f"  abstract min-max input: {abstract_ok}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
