#!/usr/bin/env python3
r"""Quotient-norm limsup theorem for augmented traces."""

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
        "--trace-recovery-json",
        default="augmented_trace_recovery_sequence_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_trace_quotient_limsup_theorem.json")
    args = parser.parse_args()

    trace = load(args.trace_recovery_json)
    trace_ok = bool(trace.get("traceRecoverySequenceClosed"))
    limsup_ok = trace_ok

    data = {
        "theoremName": "augmented trace quotient limsup theorem",
        "proofClass": "analytic proof",
        "traceRecoveryJson": args.trace_recovery_json,
        "statement": (
            "For every x=R_aug f in X_aug, finite representatives x_N=R_aug,N f_N "
            "satisfy limsup ||x_N||_{X_aug,N} <= ||x||_{X_aug}."
        ),
        "statuses": {
            "traceRecoveryInputStatus": status(
                "trace recovery input",
                trace_ok,
                "Recovery sequences recover both graph vectors and their augmented trace coordinates.",
            ),
            "traceQuotientLimsupStatus": status(
                "trace quotient limsup",
                limsup_ok,
                (
                    "Use recovered representatives f_N for x_N; the finite "
                    "quotient norm is bounded by ||f_N||_V, and ||f_N||_V -> ||f||_V "
                    "after taking near-minimizers for the limiting quotient norm."
                ),
            ),
        },
        "traceQuotientLimsupClosed": limsup_ok,
        "finiteQuotientUpperBoundClosed": limsup_ok,
        "proof": [
            "Choose f with R_aug f=x and ||f||_V within epsilon of the quotient norm.",
            "Use trace recovery to get f_N and x_N=R_aug,N f_N -> x.",
            "Bound ||x_N||_{X_aug,N} by ||f_N||_V.",
            "Pass N->infinity and epsilon->0.",
        ],
        "remainingAnalyticGap": None if limsup_ok else "Trace recovery is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace quotient limsup theorem")
    print(f"  trace recovery: {trace_ok}")
    print(f"  quotient limsup: {limsup_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
