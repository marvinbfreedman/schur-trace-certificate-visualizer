#!/usr/bin/env python3
r"""Original Weyl kernel positivity theorem for |omega|<1/2.

This theorem is the narrow input needed by the Weyl operator-family layer.  It
now imports the original Weyl positivity assembly theorem, which keeps the
branch weights, Green contraction, and exact quadratic-form identity one layer
lower in the audit graph.

The conclusion is positivity of the original coordinate Weyl kernel
K_omega for every |omega|<1/2.  It deliberately does not package the result as
the KLM quantum positive-type condition; that is handled by
``uniform_omega_weyl_klm_bridge.py`` via the hbar=1 KLM/Weyl theorem.
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


def nested_closed(data: dict, *keys: str) -> bool:
    item = data
    for key in keys:
        item = item.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--assembly-json",
        default="original_weyl_kernel_positivity_consequence_theorem.json",
    )
    parser.add_argument("--green-closure-json", default="")
    parser.add_argument("--form-transport-json", default="")
    parser.add_argument("--lift-json", default="")
    parser.add_argument("--json-out", default="original_weyl_kernel_positivity_theorem.json")
    args = parser.parse_args()

    assembly = load(args.assembly_json)
    green = load(args.green_closure_json) if args.green_closure_json else {}
    transport_path = args.form_transport_json or args.lift_json
    transport = load(transport_path) if transport_path else {}

    branch_ok = bool(
        assembly.get("uniformOmegaBranchTransportClosed")
        or assembly.get("originalWeylKernelPositivityAssemblyClosed")
    )
    green_ok = bool(
        assembly.get("greenLiftContractionClosed")
        or green.get("continuumGreenLiftClosureClosed")
        or green.get("closedOnCompletedVolterraDomain")
        or nested_closed(green, "statuses", "greenLiftContractionStatus")
    )
    transport_ok = bool(
        assembly.get("originalWeylFormTransportClosed")
        or transport.get("originalWeylFormTransportClosed")
        or transport.get("originalWeylKernelPositivityTransportClosed")
        or transport.get("defectFreePullbackClosed")
        or nested_closed(transport, "statuses", "originalFormTransportStatus")
    )
    assembly_ok = bool(
        assembly.get("originalWeylKernelPositivityConsequenceClosed")
        or assembly.get("originalWeylKernelPositivityAssemblyClosed")
        or assembly.get("originalWeylKernelPositivityClosed")
    )
    original_ok = assembly_ok and branch_ok and green_ok and transport_ok

    data = {
        "theoremName": "original Weyl kernel positivity theorem",
        "proofClass": "symbolic identity",
        "omegaRange": "|omega| < 1/2",
        "assemblyJson": args.assembly_json,
        "legacyInputs": {
            "greenClosureJson": args.green_closure_json or None,
            "formTransportJson": transport_path or None,
            "liftJson": args.lift_json or None,
        },
        "statuses": {
            "positivityAssemblyInputStatus": status(
                "original Weyl positivity assembly input",
                assembly_ok,
                "The assembly theorem packages branch weights, Green contraction, and the exact quadratic-form identity.",
            ),
            "omegaIndependentHardyMultiplierStatus": status(
                "Green multiplier contraction input",
                green_ok,
                "Imported through the original Weyl positivity assembly theorem.",
                blocker=None if green_ok else "Close continuum_green_lift_closure_theorem.py.",
            ),
            "quotientToOriginalLiftStatus": status(
                "exact original Weyl form identity input",
                transport_ok,
                "Imported through the original Weyl positivity assembly theorem.",
                blocker=None if transport_ok else "Close original_weyl_quadratic_form_identity_theorem.py.",
            ),
            "originalWeylPositivityStatus": status(
                "original Weyl kernel positivity for |omega|<1/2",
                original_ok,
                (
                    "The assembly theorem gives positivity of K_omega for every |omega|<1/2."
                ),
                blocker=None if original_ok else "Need the original Weyl positivity assembly theorem.",
            ),
        },
        "uniformOmegaBranchTransportClosed": branch_ok,
        "greenLiftContractionClosed": green_ok,
        "originalWeylFormTransportClosed": transport_ok,
        "quotientToOriginalWeylLiftClosed": transport_ok,
        "originalWeylKernelPositivityClosed": original_ok,
        "proof": [
            "Import the original Weyl positivity assembly theorem.",
            "Read off branch weights, Green contraction, and exact quadratic-form identity.",
            "Export only original coordinate Weyl kernel positivity for every |omega|<1/2.",
        ],
        "remainingAnalyticGap": None if original_ok else "One original Weyl positivity input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl kernel positivity theorem")
    print(f"  branch transport: {branch_ok}")
    print(f"  Green contraction: {green_ok}")
    print(f"  form transport: {transport_ok}")
    print(f"  original Weyl positivity: {original_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
