#!/usr/bin/env python3
r"""Trace repair descent theorem for the transported quotient X_aug."""

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
        "--descent-principle-json",
        default="abstract_bounded_form_quotient_descent_theorem.json",
    )
    parser.add_argument(
        "--null-fiber-json",
        default="augmented_repair_null_fiber_compatibility_theorem.json",
    )
    parser.add_argument(
        "--uniform-bound-json",
        default="augmented_repair_uniform_quotient_bound_theorem.json",
    )
    parser.add_argument("--quotient-json", default="")
    parser.add_argument("--json-out", default="augmented_trace_repair_descends_theorem.json")
    args = parser.parse_args()

    principle = load(args.descent_principle_json)
    null_fiber = load(args.null_fiber_json)
    uniform_bound = load(args.uniform_bound_json)

    principle_ok = bool(principle.get("abstractBoundedFormQuotientDescentClosed"))
    null_fiber_ok = bool(null_fiber.get("repairNullFiberCompatibilityClosed"))
    bound_ok = bool(uniform_bound.get("uniformRepairQuotientBoundClosed"))
    descent_ok = principle_ok and null_fiber_ok and bound_ok

    data = {
        "theoremName": "augmented trace repair descent theorem",
        "proofClass": "analytic proof",
        "descentPrincipleJson": args.descent_principle_json,
        "nullFiberJson": args.null_fiber_json,
        "uniformBoundJson": args.uniform_bound_json,
        "legacyQuotientJson": args.quotient_json or None,
        "statement": (
            "The finite augmented repair forms descend to a well-defined "
            "closed quadratic form D_aug on the transported quotient trace "
            "space X_aug."
        ),
        "statuses": {
            "abstractDescentPrincipleInputStatus": status(
                "abstract bounded-form descent input",
                principle_ok,
                (
                    "The abstract quotient descent theorem says trace-compatible "
                    "bounded forms descend uniquely to closed quotient forms."
                ),
            ),
            "nullFiberCompatibilityInputStatus": status(
                "null-fiber compatibility input",
                null_fiber_ok,
                (
                    "Finite repair null-fiber compatibility makes the repair "
                    "independent of representative in the transported trace quotient."
                ),
            ),
            "uniformQuotientBoundInputStatus": status(
                "uniform quotient bound input",
                bound_ok,
                (
                    "The interval finite repair upper bound gives a uniform "
                    "bound in the transported quotient norm."
                ),
            ),
            "traceRepairDescendsStatus": status(
                "trace repair descends to X_aug",
                descent_ok,
                (
                    "Null-fiber compatibility and the uniform quotient bound "
                    "verify the hypotheses of the abstract bounded-form "
                    "quotient descent theorem."
                ),
            ),
        },
        "traceRepairDescendsClosed": descent_ok,
        "transportedRepairFormWellDefinedClosed": descent_ok,
        "transportedRepairFormBoundedClosed": bound_ok,
        "transportedRepairRepresentativeIndependenceClosed": null_fiber_ok,
        "proof": [
            "Import finite/transported null-fiber compatibility for the repair forms.",
            "Import the uniform quotient-norm repair bound.",
            "Apply the abstract bounded-form quotient descent principle.",
            "The descended bounded form is the closed trace-side repair D_aug on X_aug.",
        ],
        "remainingAnalyticGap": None
        if descent_ok
        else "Descent principle, null-fiber compatibility, or uniform quotient bound is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented trace repair descent theorem")
    print(f"  descent principle: {principle_ok}")
    print(f"  null-fiber compatibility: {null_fiber_ok}")
    print(f"  uniform quotient bound: {bound_ok}")
    print(f"  repair descent: {descent_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
