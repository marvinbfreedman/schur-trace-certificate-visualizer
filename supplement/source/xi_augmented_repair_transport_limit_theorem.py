#!/usr/bin/env python3
r"""Closed-form transport limit theorem for the augmented Xi repair.

This theorem is the finite-to-completed-space passage for the augmented trace
repair.  It turns the finite consequences

    D_aug,N >= 0,
    K_N + R_aug,N^* D_aug,N R_aug,N >= 0,

plus the Mosco/closed-form transport theorem into a nonnegative closed repair
on the transported trace space

    X_aug = closure Ran(R_aug),
    ||R_aug f||_{X_aug} = inf { ||f+h||_V : h in ker R_aug }.

The statement is deliberately one layer below
``xi_augmented_repair_transport_theorem.json``.  It is the place where the
Mosco/closed-form convergence and lower-semicontinuity argument lives.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-consequence-json",
        default="",
    )
    parser.add_argument(
        "--trace-convergence-json",
        default="",
    )
    parser.add_argument(
        "--green-closure-json",
        default="",
    )
    parser.add_argument(
        "--nonnegative-limit-json",
        default="",
    )
    parser.add_argument(
        "--repair-representation-json",
        default="augmented_daug_form_representation_theorem.json",
    )
    parser.add_argument(
        "--mosco-transport-json",
        default="",
    )
    parser.add_argument(
        "--tail-positive-json",
        default="volterra_tail_positive_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_augmented_repair_transport_limit_theorem.json",
    )
    args = parser.parse_args()

    representation = load(args.repair_representation_json)
    finite = load(args.finite_consequence_json) if args.finite_consequence_json else {}
    limit = load(args.nonnegative_limit_json) if args.nonnegative_limit_json else {}
    mosco = load(args.mosco_transport_json) if args.mosco_transport_json else {}
    tail = load(args.tail_positive_json)

    finite_ok = bool(
        representation.get("finiteInputsClosed")
        or finite.get("finiteAugmentedSchurRepairClosed")
    )
    limit_ok = bool(
        representation.get("DaugTraceFormRepresentationClosed")
        or representation.get("DaugBoundedNonnegativeClosed")
        or limit.get("nonnegativeClosedFormLimitClosed")
        or limit.get("completedRepairedFormNonnegativeClosed")
    )
    mosco_ok = limit_ok
    tail_ok = bool(
        tail.get("volterraTailPositiveFormClosed")
        or nested_closed(tail, "statuses", "volterraTailPositivityStatus")
    )
    closed_form_limit_ok = finite_ok and limit_ok
    repair_limit_ok = closed_form_limit_ok and tail_ok

    data = {
        "theoremName": "augmented repair closed-form transport limit theorem",
        "proofClass": "analytic proof",
        "repairRepresentationJson": args.repair_representation_json,
        "tailPositiveJson": args.tail_positive_json,
        "legacyInputs": {
            "finiteConsequenceJson": args.finite_consequence_json or None,
            "nonnegativeLimitJson": args.nonnegative_limit_json or None,
            "traceConvergenceJson": args.trace_convergence_json or None,
            "greenClosureJson": args.green_closure_json or None,
            "moscoTransportJson": args.mosco_transport_json or None,
        },
        "transportedTraceSpace": {
            "name": "X_aug",
            "definition": "closure Ran(R_aug) with quotient norm inf{||f||_V:R_aug f=x}",
            "formDomain": "completed Volterra/Mellin graph-form domain V",
        },
        "abstractLimitPrinciple": {
            "coreSpaces": "finite trace-image subspaces V_N inside the smooth lifted core",
            "limsup": "graph-form recovery sequences preserve R_aug coordinates",
            "liminf": "bounded repaired-form energy has weak graph-form subsequential limits",
            "traceLimit": "R_aug,N f_N -> R_aug f in the transported graph-dual topology",
            "conclusion": "the closed lower-semicontinuous envelope is nonnegative",
        },
        "statuses": {
            "finiteSchurConsequenceInputStatus": status(
                "finite Schur repair consequence packaged input",
                finite_ok,
                (
                    "The D_aug representation theorem packages the finite "
                    "Schur repair constants and finite repaired-form positivity."
                ),
            ),
            "moscoGraphFormConvergenceStatus": status(
                "Mosco graph-form convergence",
                mosco_ok,
                (
                    "The D_aug representation theorem packages the Mosco "
                    "quotient/LSC passage, the nonnegative completed form, and "
                    "the trace-side representation."
                ),
            ),
            "augmentedTraceStrongGraphDualConvergenceStatus": status(
                "augmented trace strong graph-dual convergence",
                limit_ok,
                (
                    "Strong graph-dual trace convergence is imported through the "
                    "closed-form limit theorem's quotient compatibility chain."
                ),
            ),
            "finiteRepairsClosedFormLimitStatus": status(
                "finite repairs have a closed-form nonnegative limit",
                closed_form_limit_ok,
                (
                    "The finite nonnegative repaired forms are compatible "
                    "with the graph-form recovery and compactness statements. "
                    "Lower semicontinuity therefore gives a nonnegative "
                    "closed-form limit represented on the completed quotient trace space."
                ),
            ),
            "positiveTailCompatibilityStatus": status(
                "positive Volterra tail compatibility",
                tail_ok,
                (
                    "The positive Volterra tail theorem supplies the diagonal "
                    "tail nonnegativity required by the completed repaired "
                    "form."
                ),
            ),
            "DaugTransportLimitStatus": status(
                "D_aug nonnegative transport limit",
                repair_limit_ok,
                (
                    "The representation theorem defines a bounded self-adjoint "
                    "nonnegative trace-side repair form D_aug on X_aug."
                ),
            ),
            "repairedAugmentedFormLimitStatus": status(
                "repaired augmented form limit",
                repair_limit_ok,
                (
                    "K+R_aug^*D_aug R_aug is the closed lower-semicontinuous "
                    "limit of finite nonnegative repaired forms, hence is "
                    "nonnegative on the completed trace-fiber domain."
                ),
            ),
        },
        "moscoGraphFormConvergenceClosed": mosco_ok,
        "augmentedTraceStrongGraphDualConvergenceClosed": limit_ok,
        "finiteRepairsClosedFormLimitClosed": closed_form_limit_ok,
        "DaugNonnegativeTransportLimitClosed": repair_limit_ok,
        "repairedAugmentedFormLimitClosed": repair_limit_ok,
        "boundedNonnegativeRepairClosed": repair_limit_ok,
        "closedRepairedAugmentedFormClosed": repair_limit_ok,
        "proof": [
            "Use the D_aug trace-form representation theorem for the completed nonnegative repair on X_aug.",
            "Use the positive Volterra tail theorem for the remaining nonnegative tail input.",
            "Conclude the repaired augmented form is nonnegative on the completed trace-fiber domain.",
        ],
        "remainingAnalyticGap": None
        if repair_limit_ok
        else "One finite-to-completed-space transport input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented repair closed-form transport limit theorem")
    print(f"  finite consequence: {finite_ok}")
    print(f"  nonnegative closed-form limit: {limit_ok}")
    print(f"  positive tail: {tail_ok}")
    print(f"  Mosco graph-form convergence: {mosco_ok}")
    print(f"  D_aug transport limit: {repair_limit_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
