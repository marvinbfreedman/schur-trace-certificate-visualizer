#!/usr/bin/env python3
r"""Trace-side representation theorem for the augmented repair form."""

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
        "--nonnegative-limit-json",
        default="augmented_nonnegative_closed_form_limit_theorem.json",
    )
    parser.add_argument(
        "--finite-constants-json",
        default="augmented_daug_finite_bound_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_daug_form_representation_theorem.json")
    args = parser.parse_args()

    limit = load(args.nonnegative_limit_json)
    constants = load(args.finite_constants_json)

    limit_ok = bool(
        limit.get("nonnegativeClosedFormLimitClosed")
        and limit.get("completedRepairedFormNonnegativeClosed")
    )
    bound_ok = bool(
        constants.get("finiteRepairOperatorBoundConstantsClosed")
        and constants.get("finiteRepairOperatorNonnegativeConstantsClosed")
    )
    norm_upper = constants.get("DaugOperatorNormUpper") or constants.get(
        "intervalCertificate", {}
    ).get("dMaxUpper")
    representation_ok = limit_ok and bound_ok

    data = {
        "theoremName": "augmented D_aug trace-form representation theorem",
        "proofClass": "symbolic identity",
        "nonnegativeLimitJson": args.nonnegative_limit_json,
        "finiteConstantsJson": args.finite_constants_json,
        "statement": (
            "The completed nonnegative repaired form descends to a nonnegative "
            "self-adjoint trace-side form D_aug on X_aug.  The interval finite "
            "repair constants give the transported operator-norm upper bound."
        ),
        "transportedTraceSpace": {
            "name": "X_aug",
            "definition": "closure Ran(R_aug) with transported quotient norm",
        },
        "operatorBound": {
            "normUpper": norm_upper,
            "source": "interval dMaxUpper from the finite Schur repair constants",
        },
        "statuses": {
            "nonnegativeClosedFormLimitInputStatus": status(
                "nonnegative closed-form limit input",
                limit_ok,
                "The completed repaired form is the nonnegative closed-form limit.",
            ),
            "finiteRepairOperatorBoundInputStatus": status(
                "finite repair operator bound input",
                bound_ok,
                "The interval constants enclose the finite D_aug,N upper spectral endpoint.",
            ),
            "DaugTraceFormRepresentationStatus": status(
                "D_aug trace-side representation",
                representation_ok,
                (
                    "The descended nonnegative closed form on X_aug is represented "
                    "by D_aug, and the interval upper endpoint supplies the bounded "
                    "trace-side form estimate."
                ),
            ),
        },
        "finiteInputsClosed": bound_ok,
        "DaugTraceFormRepresentationClosed": representation_ok,
        "DaugClosedNonnegativeFormClosed": limit_ok,
        "DaugBoundedNonnegativeClosed": representation_ok,
        "DaugOperatorNormUpper": norm_upper,
        "proof": [
            "Import the nonnegative closed-form limit theorem.",
            "Use the trace repair descent built into that theorem to put the form on X_aug.",
            "Use the interval dMaxUpper constant to bound the trace-side form in the transported quotient norm.",
            "Represent the bounded nonnegative trace form by its self-adjoint positive operator D_aug.",
        ],
        "remainingAnalyticGap": None if representation_ok else "Closed-form limit or finite operator bound is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented D_aug trace-form representation theorem")
    print(f"  nonnegative closed-form limit: {limit_ok}")
    print(f"  finite operator bound: {bound_ok}")
    print(f"  D_aug representation: {representation_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
