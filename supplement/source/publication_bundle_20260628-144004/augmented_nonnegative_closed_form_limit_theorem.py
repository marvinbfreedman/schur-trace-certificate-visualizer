#!/usr/bin/env python3
r"""Nonnegative closed-form limit theorem for the augmented repair."""

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
        "--pointwise-nonnegative-json",
        default="augmented_pointwise_repaired_form_nonnegative_certificate.json",
    )
    parser.add_argument("--finite-sequence-json", default="")
    parser.add_argument(
        "--closed-form-lsc-json",
        default="closed_form_lsc_transport_consequence_theorem.json",
    )
    parser.add_argument(
        "--cone-principle-json",
        default="closed_lsc_nonnegative_cone_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="augmented_nonnegative_closed_form_limit_theorem.json")
    args = parser.parse_args()

    pointwise = load(args.pointwise_nonnegative_json)
    lsc = load(args.closed_form_lsc_json)
    cone = load(args.cone_principle_json)

    pointwise_ok = bool(pointwise.get("pointwiseRepairedFormsNonnegativeClosed"))
    envelope_ok = bool(
        lsc.get("closedLowerEnvelopeIdentifiedClosed")
        and (lsc.get("closedFormTransportClosed") or lsc.get("lowerSemicontinuityClosed"))
    )
    cone_ok = bool(cone.get("closedLscNonnegativeConePrincipleClosed"))
    limit_ok = pointwise_ok and envelope_ok and cone_ok

    data = {
        "theoremName": "augmented nonnegative closed-form limit theorem",
        "proofClass": "analytic proof",
        "pointwiseNonnegativeJson": args.pointwise_nonnegative_json,
        "closedFormLscJson": args.closed_form_lsc_json,
        "conePrincipleJson": args.cone_principle_json,
        "legacyFiniteSequenceJson": args.finite_sequence_json or None,
        "statement": (
            "Pointwise nonnegative augmented repaired core forms, together "
            "with an identified closed lower-semicontinuous transported "
            "envelope, produce a nonnegative completed repaired form by the "
            "abstract closed-LSC cone principle."
        ),
        "statuses": {
            "pointwiseNonnegativeCoreInputStatus": status(
                "pointwise nonnegative core input",
                pointwise_ok,
                (
                    "The pointwise certificate proves q_N^rep>=0 for each "
                    "fixed augmented trace-core index N and makes no continuum "
                    "promotion claim."
                ),
            ),
            "closedFormLscInputStatus": status(
                "closed lower-envelope transport input",
                envelope_ok,
                (
                    "The closed-form LSC transport theorem identifies the "
                    "completed lower envelope and transports it to X_aug."
                ),
            ),
            "closedLscConePrincipleInputStatus": status(
                "closed LSC cone principle input",
                cone_ok,
                (
                    "The abstract cone principle says a closed LSC lower "
                    "envelope of nonnegative core forms is nonnegative."
                ),
            ),
            "nonnegativeClosedFormLimitStatus": status(
                "nonnegative closed-form limit",
                limit_ok,
                (
                    "The pointwise core certificate supplies nonnegative core "
                    "forms; the transport theorem identifies the closed lower "
                    "envelope; the abstract cone principle keeps that envelope "
                    "inside the nonnegative cone."
                ),
            ),
        },
        "nonnegativeClosedFormLimitClosed": limit_ok,
        "completedRepairedFormNonnegativeClosed": limit_ok,
        "proof": [
            "Import pointwise nonnegativity on each fixed augmented trace core.",
            "Import the closed-form LSC transport theorem to identify the completed lower envelope on X_aug.",
            "Apply the abstract closed-LSC nonnegative cone principle to the transported envelope.",
            "Conclude the completed repaired form is nonnegative.",
        ],
        "remainingAnalyticGap": None
        if limit_ok
        else "Pointwise core nonnegativity, closed-envelope transport, or the abstract cone principle is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented nonnegative closed-form limit theorem")
    print(f"  pointwise nonnegative core: {pointwise_ok}")
    print(f"  closed-envelope transport: {envelope_ok}")
    print(f"  closed-LSC cone principle: {cone_ok}")
    print(f"  nonnegative limit: {limit_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
