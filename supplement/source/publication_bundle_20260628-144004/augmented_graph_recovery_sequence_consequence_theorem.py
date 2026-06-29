#!/usr/bin/env python3
r"""Narrow consequence of augmented graph recovery."""

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
        "--graph-recovery-json",
        default="augmented_graph_recovery_sequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_graph_recovery_sequence_consequence_theorem.json",
    )
    args = parser.parse_args()

    graph = load(args.graph_recovery_json)
    graph_ok = bool(graph.get("graphRecoverySequenceClosed"))
    finite_core_ok = bool(graph.get("finiteCoreRecoveryClosed", graph_ok))
    ok = graph_ok and finite_core_ok

    data = {
        "theoremName": "augmented graph recovery sequence consequence theorem",
        "proofClass": "symbolic identity",
        "graphRecoveryJson": args.graph_recovery_json,
        "statement": (
            "The augmented graph recovery theorem supplies the graph-norm "
            "recovery sequence consequence needed by augmented trace recovery."
        ),
        "statuses": {
            "graphRecoverySequenceStatus": status(
                "graph recovery sequence",
                graph_ok,
                "Imported from the augmented graph recovery sequence theorem.",
            ),
            "finiteCoreRecoveryStatus": status(
                "finite core recovery",
                finite_core_ok,
                "Imported from the augmented graph recovery sequence theorem.",
            ),
            "graphRecoveryConsequenceStatus": status(
                "graph recovery consequence",
                ok,
                "Finite core functions recover completed graph-domain vectors.",
            ),
        },
        "graphRecoverySequenceClosed": graph_ok,
        "finiteCoreRecoveryClosed": finite_core_ok,
        "graphRecoverySequenceConsequenceClosed": ok,
        "proof": [
            "Import the augmented graph recovery sequence theorem.",
            "Extract the graph-norm recovery conclusion.",
            "Expose only that consequence to the augmented trace recovery theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Graph recovery input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented graph recovery sequence consequence theorem")
    print(f"  graph recovery: {graph_ok}")
    print(f"  finite core recovery: {finite_core_ok}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
