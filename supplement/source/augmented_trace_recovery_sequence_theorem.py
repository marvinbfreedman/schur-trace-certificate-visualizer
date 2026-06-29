#!/usr/bin/env python3
r"""Trace recovery sequence theorem for augmented Mosco transport."""

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
        default="augmented_graph_recovery_sequence_consequence_theorem.json",
    )
    parser.add_argument(
        "--trace-convergence-json",
        default="xi_augmented_trace_convergence_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_trace_recovery_sequence_theorem.json")
    args = parser.parse_args()

    graph = load(args.graph_recovery_json)
    trace = load(args.trace_convergence_json)

    graph_ok = bool(graph.get("graphRecoverySequenceClosed"))
    trace_ok = bool(
        trace.get("closedTraceConvergenceClosed")
        and trace.get("augmentedTraceGraphConvergenceClosed")
    )
    recovery_ok = graph_ok and trace_ok

    data = {
        "theoremName": "augmented trace recovery sequence theorem",
        "proofClass": "analytic proof",
        "graphRecoveryJson": args.graph_recovery_json,
        "traceConvergenceJson": args.trace_convergence_json,
        "statement": (
            "The graph-norm recovery sequence can be chosen so that "
            "R_aug,N f_N -> R_aug f in the transported graph-dual trace topology."
        ),
        "statuses": {
            "graphRecoveryInputStatus": status(
                "graph recovery input",
                graph_ok,
                "Finite core functions recover f in the completed graph norm.",
            ),
            "traceConvergenceInputStatus": status(
                "augmented trace convergence input",
                trace_ok,
                "R_aug,N converges strongly to R_aug on graph-norm bounded recovery sequences.",
            ),
            "traceRecoverySequenceStatus": status(
                "trace recovery sequence",
                recovery_ok,
                (
                    "Graph recovery gives f_N -> f in V; strong graph-dual "
                    "trace convergence then gives R_aug,N f_N -> R_aug f."
                ),
            ),
        },
        "traceRecoverySequenceClosed": recovery_ok,
        "graphAndTraceRecoveryClosed": recovery_ok,
        "proof": [
            "Take the graph recovery sequence from the graph recovery theorem.",
            "Apply strong graph-dual convergence of the augmented trace rows.",
            "Use the triangle inequality in the transported trace norm.",
            "Conclude simultaneous recovery of f and R_aug f.",
        ],
        "remainingAnalyticGap": None if recovery_ok else "Graph recovery or trace convergence is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace recovery sequence theorem")
    print(f"  graph recovery: {graph_ok}")
    print(f"  trace convergence: {trace_ok}")
    print(f"  trace recovery: {recovery_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
