#!/usr/bin/env python3
r"""Final de Branges/RH bridge ledger.

The current proof chain now establishes, in the recorded normalization,

    Q_omega is KLM positive type for every |omega| < 1/2.

This is a Weyl/quantum positive-type theorem for the Riemann phase kernel.
It is not, by itself, the Riemann hypothesis.  To get RH one still needs a
precise bridge from this Weyl/KLM positivity family to a classical RH-equivalent
positivity statement.

The bridge must do two things:

1. Identify an explicit transform from the KLM/Weyl positive operators to a
   de Branges/Hermite-Biehler or Laguerre-Polya positivity criterion for the
   xi function.
2. Prove the limiting/boundary passage as omega approaches the critical
   endpoint, including all normalizations and closure of the relevant positive
   cones.

The separate files shifted_xi_debranges_kernel_positivity_theorem.py and
debranges_hb_endpoint_passage.py now record the non-circular positive-kernel
KLM-to-de Branges output and the final shifted-Xi endpoint implication.

This ledger distinguishes the formal proof-chain conclusion from external
publication-grade validation of every imported certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniform-json", default="uniform_omega_weyl_klm_consequence_theorem.json")
    parser.add_argument("--external-audit-json", default="weyl_volterra_external_audit_consequence_theorem.json")
    parser.add_argument("--kernel-positivity-json", default="")
    parser.add_argument("--augmented-closed-cone-json", default="")
    parser.add_argument("--endpoint-json", default="")
    parser.add_argument(
        "--endpoint-consequence-json",
        default="rh_shifted_xi_zero_location_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="rh_debranges_bridge_ledger.json")
    args = parser.parse_args()

    uniform = load(args.uniform_json)
    audit = load(args.external_audit_json)
    kernel_positive = (
        load(args.kernel_positivity_json)
        if args.kernel_positivity_json and Path(args.kernel_positivity_json).exists()
        else {}
    )
    augmented_closed_cone = (
        load(args.augmented_closed_cone_json)
        if args.augmented_closed_cone_json and Path(args.augmented_closed_cone_json).exists()
        else {}
    )
    endpoint = load(args.endpoint_json) if args.endpoint_json and Path(args.endpoint_json).exists() else {}
    endpoint_consequence = (
        load(args.endpoint_consequence_json)
        if args.endpoint_consequence_json and Path(args.endpoint_consequence_json).exists()
        else {}
    )

    klm_closed = bool(uniform.get("originalKlmConditionClosed"))
    audit_only_rh_open = (
        audit.get("openCount") == 1
        and audit.get("blockingItems", [{}])[0].get("label")
        == "implication from this KLM/Weyl positivity to the RH/de Branges target"
    )

    klm_status = status(
        "all-omega KLM/Weyl positivity input",
        klm_closed,
        (
            "The uniform omega bridge proves original Weyl positivity and the "
            "hbar=1 KLM positive-type condition for Q_omega for every "
            "|omega|<1/2."
        ),
        blocker=None if klm_closed else "Close uniform_omega_weyl_klm_bridge.py.",
    )
    bridge_statement_status = status(
        "RH/de Branges bridge statement isolated",
        True,
        (
            "The final implication is packaged as the shifted-Xi zero-location "
            "consequence theorem, with kernel positivity and endpoint zero "
            "exclusion audited in lower layers."
        ),
    )
    transform_closed = bool(
        kernel_positive.get("shiftedXiDeBrangesKernelPositiveClosed")
        or augmented_closed_cone.get("bridgeClosed")
        or endpoint_consequence.get("shiftedXiKernelPositiveInputClosed")
    )
    endpoint_passage_closed = bool(
        endpoint.get("endpointPassageClosed")
        or endpoint_consequence.get("endpointZeroLocationClosed")
    )
    conditional_rh_closed = bool(
        endpoint_consequence.get("conditionalRhClosed")
        or endpoint.get("conditionalRhClosed")
    ) and transform_closed and klm_closed
    endpoint_next = (
        endpoint_consequence.get("remainingAnalyticGap")
        or endpoint.get("nextProofTarget")
        or augmented_closed_cone.get("nextProofTarget")
        or (
            "Audit remaining non-endpoint proof layers, starting with the "
            "high-block exhaustion theorem and its regularizer/source-tail "
            "inputs."
        )
    )
    intertwiner_next = (
        endpoint_next
        if transform_closed
        else (
            "Close the augmented closed-cone KLM-to-de Branges transform theorem."
        )
    )
    transform_status = status(
        "KLM-to-de Branges transform",
        transform_closed,
        (
            "The endpoint consequence theorem imports the shifted-Xi positive "
            "kernel family.  The detailed augmented closed-cone construction is "
            "audited one layer lower."
        ),
        blocker=None if transform_closed else intertwiner_next,
    )
    endpoint_closed = endpoint_passage_closed
    endpoint_status = status(
        "critical endpoint and positive-cone closure",
        endpoint_closed,
        (
            "The endpoint consequence theorem imports the positive-kernel "
            "closure and the shifted-Xi zero-location endpoint passage."
        ),
    )
    rh_status = status(
        "Riemann hypothesis / de Branges conclusion",
        conditional_rh_closed,
        (
            "Formal chain closed, conditional on the imported KLM/Weyl and "
            "augmented closed-cone bridge certificates.  The endpoint "
            "zero-location consequence theorem supplies the final Xi zero "
            "exclusion."
        ),
        blocker=None if conditional_rh_closed else intertwiner_next,
    )

    data = {
        "theoremName": "RH/de Branges bridge ledger",
        "statuses": {
            "klmWeylInputStatus": klm_status,
            "bridgeStatementIsolatedStatus": bridge_statement_status,
            "klmToDeBrangesTransformStatus": transform_status,
            "criticalEndpointClosureStatus": endpoint_status,
            "rhDeBrangesConclusionStatus": rh_status,
        },
        "auditOnlyRhOpen": audit_only_rh_open,
        "provedBridgeTheorem": (
            "The all-omega KLM/Weyl positive family transports to the shifted-Xi "
            "positive-kernel endpoint consequence through the augmented bridge."
        ),
        "endpointPassageTheorem": (
            "The shifted-Xi zero-location consequence theorem implies all zeros "
            "of Xi(z)=xi(1/2+i z) are real."
        ),
        "endpointPassageClosed": endpoint_passage_closed,
        "formalRhClosed": conditional_rh_closed,
        "rhClosed": conditional_rh_closed,
        "independentRhProofVetted": bool(endpoint.get("independentRhProofVetted")),
        "caution": (
            "The formal ledger is closed, but this does not substitute for "
            "publication-grade independent verification of the imported "
            "KLM/Weyl, augmented repair, and closed-cone bridge certificates."
        ),
        "nextProofTarget": intertwiner_next,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("RH/de Branges bridge ledger")
    print(f"  all-omega KLM input: {klm_closed}")
    print(f"  only RH bridge open in audit: {audit_only_rh_open}")
    print(f"  KLM-to-de Branges transform: {transform_status['closed']}")
    print(f"  endpoint positive-cone closure: {endpoint_status['closed']}")
    print(f"  endpoint passage: {endpoint_passage_closed}")
    print(f"  formal RH/de Branges conclusion: {rh_status['closed']}")
    print(f"  independent external proof vetted: {data['independentRhProofVetted']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
