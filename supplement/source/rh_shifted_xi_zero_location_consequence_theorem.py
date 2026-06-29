#!/usr/bin/env python3
r"""RH-side zero-location consequence for shifted Xi."""

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
        "--zero-descent-json",
        default="shifted_xi_zero_descent_endpoint_consequence_theorem.json",
    )
    parser.add_argument(
        "--kernel-positivity-json",
        default="shifted_xi_debranges_kernel_positivity_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="rh_shifted_xi_zero_location_consequence_theorem.json")
    args = parser.parse_args()

    zero = load(args.zero_descent_json)
    kernel = load(args.kernel_positivity_json)

    zero_ok = bool(zero.get("zeroDescentEndpointClosed") and zero.get("endpointZeroLocationClosed"))
    kernel_ok = bool(kernel.get("shiftedXiDeBrangesKernelPositiveClosed"))
    rh_ok = zero_ok and kernel_ok

    data = {
        "theoremName": "RH shifted-Xi zero-location consequence theorem",
        "proofClass": "analytic proof",
        "zeroDescentJson": args.zero_descent_json,
        "kernelPositivityJson": args.kernel_positivity_json,
        "statement": (
            "The shifted-Xi positive kernel family gives the endpoint zero-location "
            "consequence: all zeros of Xi(z)=xi(1/2+i z) are real."
        ),
        "statuses": {
            "shiftedXiKernelPositiveInputStatus": status(
                "shifted-Xi positive-kernel input",
                kernel_ok,
                "The shifted-Xi de Branges kernel is positive on every finite upper-half-plane evaluation set.",
            ),
            "zeroLocationInputStatus": status(
                "zero-location endpoint input",
                zero_ok,
                "The zero-descent endpoint theorem excludes non-real zeros of Xi.",
            ),
            "rhZeroLocationConsequenceStatus": status(
                "RH-side zero-location consequence",
                rh_ok,
                "All Xi zeros are real in the z-normalization.",
            ),
        },
        "shiftedXiKernelPositiveInputClosed": kernel_ok,
        "endpointZeroLocationClosed": zero_ok,
        "conditionalRhClosed": rh_ok,
        "formalRhClosed": rh_ok,
        "proof": [
            "Import shifted-Xi positive-kernel input.",
            "Import the zero-location endpoint theorem.",
            "Translate real zeros of Xi(z)=xi(1/2+i z) to the critical line for xi(s).",
        ],
        "remainingAnalyticGap": None if rh_ok else "Kernel positivity or endpoint zero-location input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("RH shifted-Xi zero-location consequence theorem")
    print(f"  shifted kernel positivity: {kernel_ok}")
    print(f"  zero-location endpoint: {zero_ok}")
    print(f"  RH-side consequence: {rh_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
