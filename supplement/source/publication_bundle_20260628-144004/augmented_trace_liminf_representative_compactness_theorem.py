#!/usr/bin/env python3
r"""Representative compactness for the augmented trace quotient liminf step.

This theorem supplies the model-specific hypotheses required by the abstract
quotient-liminf principle.  It does not itself perform the quotient-norm
inequality.
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
    parser.add_argument("--mosco-liminf-json", default="augmented_mosco_liminf_theorem.json")
    parser.add_argument(
        "--trace-limit-json",
        default="xi_augmented_trace_limit_identification_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_trace_liminf_representative_compactness_theorem.json",
    )
    args = parser.parse_args()

    liminf = load(args.mosco_liminf_json)
    trace = load(args.trace_limit_json)

    mosco_ok = bool(
        liminf.get("moscoLiminfClosed")
        and liminf.get("weakGraphCompactnessClosed")
        and liminf.get("graphFormLowerSemicontinuityClosed")
    )
    trace_ok = bool(
        trace.get("traceLimitIdentificationClosed")
        or (
            trace.get("closedTraceConvergenceClosed")
            and trace.get("augmentedTraceGraphConvergenceClosed")
        )
    )
    representative_ok = mosco_ok and trace_ok

    data = {
        "theoremName": "augmented trace liminf representative compactness theorem",
        "proofClass": "analytic proof",
        "moscoLiminfJson": args.mosco_liminf_json,
        "traceLimitJson": args.trace_limit_json,
        "statement": (
            "Near-minimizing finite representatives for bounded augmented "
            "trace quotient coordinates have weak graph subsequential limits, "
            "and the augmented trace of the limit is the limiting trace coordinate."
        ),
        "statuses": {
            "weakRepresentativeCompactnessInputStatus": status(
                "weak representative compactness input",
                mosco_ok,
                (
                    "The Mosco liminf theorem supplies weak graph compactness "
                    "and graph-norm lower semicontinuity for bounded "
                    "representative sequences."
                ),
            ),
            "traceLimitIdentificationInputStatus": status(
                "trace limit identification input",
                trace_ok,
                (
                    "The augmented trace convergence theorem identifies the "
                    "trace of the weak graph limit with the limiting augmented "
                    "trace coordinate."
                ),
            ),
            "representativeCompactnessAndTraceIdentificationStatus": status(
                "representative compactness and trace identification",
                representative_ok,
                (
                    "Near-minimizers are bounded in the graph norm, hence have "
                    "weak graph limit points, and the closed augmented trace "
                    "map identifies the limit trace."
                ),
            ),
        },
        "representativeCompactnessAndTraceIdentificationClosed": representative_ok,
        "weakRepresentativeCompactnessClosed": mosco_ok,
        "traceLimitIdentificationClosed": trace_ok,
        "proof": [
            "Choose finite representatives within epsilon_N of their finite quotient norms.",
            "Bounded quotient norms give bounded representative graph norms.",
            "Apply Mosco liminf weak compactness to extract a weak graph limit.",
            "Apply closed augmented trace convergence to identify the limit coordinate.",
            "Retain graph-norm lower semicontinuity for the abstract quotient-liminf principle.",
        ],
        "remainingAnalyticGap": None
        if representative_ok
        else "Mosco liminf compactness or augmented trace convergence is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace liminf representative compactness theorem")
    print(f"  Mosco liminf compactness: {mosco_ok}")
    print(f"  trace limit identification: {trace_ok}")
    print(f"  representative compactness: {representative_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
