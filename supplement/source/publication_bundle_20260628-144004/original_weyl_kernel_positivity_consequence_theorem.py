#!/usr/bin/env python3
r"""Narrow consequence of original Weyl positivity assembly.

This file exports only the sign consequence needed by the Weyl operator-family
layer:

    K_omega >= 0 for every |omega| < 1/2.

The branch-weight, Green-contraction, and quadratic-form identity inputs remain
in ``original_weyl_kernel_positivity_assembly_theorem.json``.
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
    parser.add_argument(
        "--assembly-json",
        default="original_weyl_kernel_positivity_assembly_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="original_weyl_kernel_positivity_consequence_theorem.json",
    )
    args = parser.parse_args()

    assembly = load(args.assembly_json)
    assembly_ok = bool(
        assembly.get("originalWeylKernelPositivityAssemblyClosed")
        or assembly.get("originalWeylKernelPositivityClosed")
    )
    branch_ok = bool(assembly.get("uniformOmegaBranchTransportClosed"))
    green_ok = bool(assembly.get("greenLiftContractionClosed"))
    identity_ok = bool(assembly.get("originalWeylFormTransportClosed"))
    ok = assembly_ok and branch_ok and green_ok and identity_ok

    data = {
        "theoremName": "original Weyl kernel positivity consequence theorem",
        "proofClass": "symbolic identity",
        "assemblySource": "original Weyl kernel positivity assembly theorem",
        "omegaRange": "|omega| < 1/2",
        "statuses": {
            "assemblyInputStatus": status(
                "original Weyl positivity assembly input",
                assembly_ok,
                "The assembly theorem proves original Weyl kernel positivity.",
            ),
            "branchWeightSummaryStatus": status(
                "branch-weight summary",
                branch_ok,
                "The assembly theorem includes the branch-weight input.",
            ),
            "greenContractionSummaryStatus": status(
                "Green contraction summary",
                green_ok,
                "The assembly theorem includes the Green contraction input.",
            ),
            "quadraticIdentitySummaryStatus": status(
                "quadratic identity summary",
                identity_ok,
                "The assembly theorem includes the exact quadratic-form identity input.",
            ),
            "originalWeylKernelPositivityConsequenceStatus": status(
                "original Weyl kernel positivity consequence",
                ok,
                "The only exported consequence is K_omega>=0 for every |omega|<1/2.",
            ),
        },
        "uniformOmegaBranchTransportClosed": branch_ok,
        "greenLiftContractionClosed": green_ok,
        "originalWeylFormTransportClosed": identity_ok,
        "originalWeylKernelPositivityConsequenceClosed": ok,
        "originalWeylKernelPositivityClosed": ok,
        "proof": [
            "Import the original Weyl kernel positivity assembly theorem.",
            "Export only the positivity conclusion K_omega>=0 for |omega|<1/2.",
        ],
        "notExportedHere": [
            "branch-weight proof",
            "Green contraction proof",
            "quadratic-form identity proof",
        ],
        "remainingAnalyticGap": None if ok else "Original Weyl positivity assembly input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl kernel positivity consequence theorem")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
