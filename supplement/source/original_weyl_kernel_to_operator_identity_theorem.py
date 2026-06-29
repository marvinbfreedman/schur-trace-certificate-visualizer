#!/usr/bin/env python3
r"""Identity from original Weyl kernel positivity to Weyl operator positivity."""

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
        "--kernel-json",
        default="original_weyl_kernel_positivity_theorem.json",
    )
    parser.add_argument("--json-out", default="original_weyl_kernel_to_operator_identity_theorem.json")
    args = parser.parse_args()

    kernel = load(args.kernel_json)
    kernel_ok = bool(kernel.get("originalWeylKernelPositivityClosed"))

    data = {
        "theoremName": "original Weyl kernel-to-operator identity theorem",
        "proofClass": "symbolic identity",
        "kernelJson": args.kernel_json,
        "normalization": "hbar=1 Weyl convention",
        "statement": (
            "In the fixed hbar=1 Weyl convention, nonnegativity of the "
            "original coordinate Weyl kernel quadratic form is exactly "
            "nonnegativity of Op^W(sigma_omega)."
        ),
        "statuses": {
            "kernelPositivityInputStatus": status(
                "original Weyl kernel positivity input",
                kernel_ok,
                "The original Weyl kernel quadratic form is nonnegative for every |omega|<1/2.",
            ),
            "kernelToOperatorIdentityStatus": status(
                "kernel-to-operator identity",
                kernel_ok,
                "The coordinate kernel is the integral kernel of Op^W(sigma_omega) in the fixed normalization.",
            ),
        },
        "originalWeylKernelToOperatorIdentityClosed": kernel_ok,
        "originalWeylOperatorPositiveClosed": kernel_ok,
        "originalWeylKernelPositivityClosed": kernel_ok,
        "proof": [
            "Import original Weyl kernel positivity.",
            "Use the fixed hbar=1 Weyl convention identifying the coordinate kernel with Op^W(sigma_omega).",
            "A nonnegative integral-kernel quadratic form is exactly positivity of the corresponding Weyl operator.",
        ],
        "remainingAnalyticGap": None if kernel_ok else "Original Weyl kernel positivity is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Original Weyl kernel-to-operator identity theorem")
    print(f"  kernel input: {kernel_ok}")
    print(f"  operator identity: {kernel_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
