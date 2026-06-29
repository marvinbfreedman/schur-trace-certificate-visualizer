#!/usr/bin/env python3
r"""Mosco limsup recovery theorem for augmented trace transport."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--graph-recovery-json",
        default="augmented_graph_recovery_sequence_theorem.json",
    )
    parser.add_argument(
        "--trace-recovery-json",
        default="augmented_trace_recovery_sequence_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_mosco_limsup_theorem.json")
    args = parser.parse_args()

    graph = load(args.graph_recovery_json)
    trace = load(args.trace_recovery_json)

    graph_ok = bool(graph.get("graphRecoverySequenceClosed"))
    trace_ok = bool(trace.get("traceRecoverySequenceClosed"))
    limsup_ok = graph_ok and trace_ok

    data = {
        "theoremName": "augmented Mosco limsup theorem",
        "proofClass": "analytic proof",
        "graphRecoveryJson": args.graph_recovery_json,
        "traceRecoveryJson": args.trace_recovery_json,
        "statement": (
            "For every f in V there are finite core functions f_N with "
            "f_N -> f in the augmented graph norm and R_aug,N f_N -> R_aug f "
            "in the transported graph-dual trace topology."
        ),
        "statuses": {
            "graphRecoveryInputStatus": status(
                "graph recovery input",
                graph_ok,
                "Finite smooth core functions recover f in the completed graph norm.",
            ),
            "traceRecoveryInputStatus": status(
                "trace recovery input",
                trace_ok,
                "The same recovery sequence recovers R_aug f in graph-dual trace topology.",
            ),
            "moscoLimsupStatus": status(
                "Mosco limsup recovery",
                limsup_ok,
                (
                    "Graph recovery supplies f_N -> f; trace recovery supplies "
                    "R_aug,N f_N -> R_aug f.  Together these are exactly the "
                    "Mosco limsup recovery statement."
                ),
            ),
        },
        "moscoLimsupClosed": limsup_ok,
        "graphRecoverySequenceClosed": graph_ok,
        "traceRecoverySequenceClosed": trace_ok,
        "proof": [
            "Import the graph-norm recovery sequence theorem.",
            "Import the trace recovery sequence theorem.",
            "Combine the two recovered convergences in the product graph/trace topology.",
        ],
        "remainingAnalyticGap": None if limsup_ok else "Graph recovery or trace recovery is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented Mosco limsup theorem")
    print(f"  graph recovery: {graph_ok}")
    print(f"  trace recovery: {trace_ok}")
    print(f"  limsup: {limsup_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
