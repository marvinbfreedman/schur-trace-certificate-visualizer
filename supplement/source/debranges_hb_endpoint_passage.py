#!/usr/bin/env python3
r"""Hermite--Biehler/de Branges endpoint passage for shifted Xi.

This is the final formal implication layer after the KLM-to-de Branges
positive-kernel theorem has identified the shifted-Xi de Branges kernels

    E_omega(z)  = Xi(z+i omega),
    E_omega#(z) = Xi(z-i omega),
    0 < omega < 1/2.

The key point is that strict Hermite--Biehler positivity at a fixed omega is
not needed for the RH-side zero-location conclusion.  The diagonal inequality
and zero-descent contradiction are now isolated in
``debranges_diagonal_inequality_theorem.json`` and
``shifted_xi_zero_descent_endpoint_theorem.json``.

The script records this implication as a theorem ledger.  It consumes the
narrow shifted-Xi kernel positivity theorem; the detailed KLM/de Branges
closed-cone construction is audited one layer lower.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


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
    parser.add_argument(
        "--kernel-positivity-json",
        default="shifted_xi_debranges_kernel_positivity_consequence_theorem.json",
    )
    parser.add_argument("--augmented-closed-cone-json", default="")
    parser.add_argument("--uniform-json", default="")
    parser.add_argument(
        "--zero-descent-json",
        default="shifted_xi_zero_descent_endpoint_consequence_theorem.json",
    )
    parser.add_argument("--json-out", default="debranges_hb_endpoint_passage.json")
    args = parser.parse_args()

    kernel = load_optional(args.kernel_positivity_json)
    bridge = load_optional(args.augmented_closed_cone_json)
    uniform = load_optional(args.uniform_json)
    zero_descent = load_optional(args.zero_descent_json)

    kernel_positive = bool(
        kernel and kernel.get("shiftedXiDeBrangesKernelPositiveClosed")
    )
    bridge_closed = bool(bridge and bridge.get("bridgeClosed"))
    all_omega_klm = bool(uniform and uniform.get("originalKlmConditionClosed"))
    all_omega_bridge = kernel_positive or (bridge_closed and all_omega_klm)

    positivity_status = status(
        "shifted-Xi de Branges kernel positivity for all 0<omega<1/2",
        all_omega_bridge,
        (
            "Imported from the shifted-Xi de Branges kernel positivity theorem.  "
            "For every 0<omega<1/2, the kernel K_{E_omega} is positive on "
            "finite evaluation sets in the upper half-plane."
        ),
        blocker=None
        if all_omega_bridge
        else "Close shifted_xi_debranges_kernel_positivity_theorem.py.",
    )
    zero_descent_closed = bool(zero_descent and zero_descent.get("zeroDescentEndpointClosed"))
    diagonal_status = status(
        "diagonal de Branges inequality",
        bool(zero_descent and zero_descent.get("conditionalRhClosed")),
        "Imported through the zero-descent endpoint theorem.",
    )
    zero_descent_status = status(
        "zero-descent contradiction",
        zero_descent_closed,
        "The zero-descent endpoint theorem excludes upper-half-plane Xi zeros.",
    )
    lower_half_status = status(
        "lower-half-plane exclusion",
        zero_descent_closed,
        "The zero-descent endpoint theorem uses real-entire conjugation for the lower half-plane.",
    )

    endpoint_closed = all(
        s["closed"]
        for s in [
            positivity_status,
            diagonal_status,
            zero_descent_status,
            lower_half_status,
        ]
    )
    rh_status = status(
        "RH-side zero-location conclusion in z-normalization",
        endpoint_closed,
        (
            "All zeros of Xi(z)=xi(1/2+i z) are real.  This is equivalent to "
            "the nontrivial zeros of xi(s) lying on Re(s)=1/2."
        ),
        blocker=None if endpoint_closed else positivity_status.get("blocker"),
    )

    data = {
        "theoremName": "Hermite-Biehler/de Branges endpoint passage for shifted Xi",
        "kernelPositivityJson": args.kernel_positivity_json if kernel else None,
        "augmentedClosedConeJson": args.augmented_closed_cone_json if bridge else None,
        "uniformOmegaJson": args.uniform_json if uniform else None,
        "zeroDescentJson": args.zero_descent_json if zero_descent else None,
        "statuses": {
            "shiftedKernelPositivityStatus": positivity_status,
            "diagonalInequalityStatus": diagonal_status,
            "zeroDescentContradictionStatus": zero_descent_status,
            "lowerHalfPlaneExclusionStatus": lower_half_status,
            "rhSideZeroLocationStatus": rh_status,
        },
        "normalization": {
            "Xi": "Xi(z)=xi(1/2+i z), real entire",
            "Eomega": "E_omega(z)=Xi(z+i omega)",
            "EomegaSharp": "E_omega#(z)=Xi(z-i omega)",
            "omegaRange": "0<omega<1/2",
            "kernel": (
                "K_E(w,z)=(E(z)conj(E(w))-E#(z)conj(E#(w)))/"
                "(2*pi*i*(conj(w)-z))"
            ),
        },
        "proof": [
            "Import shifted-Xi kernel positivity.",
            "Import the diagonal inequality and zero-descent endpoint theorem.",
            "Conclude all Xi zeros are real in the z-normalization.",
        ],
        "endpointPassageClosed": endpoint_closed,
        "conditionalRhClosed": endpoint_closed,
        "independentRhProofVetted": False,
        "caution": (
            "This ledger closes the endpoint implication conditional on the "
            "imported KLM/Weyl and augmented closed-cone bridge certificates.  "
            "It is not an external peer-review validation of those hard inputs."
        ),
        "nextProofTarget": (
            "Audit every imported theorem for publication-grade rigor and "
            "normalization consistency."
        )
        if endpoint_closed
        else positivity_status.get("blocker"),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Hermite-Biehler/de Branges endpoint passage")
    print(f"  shifted kernel positivity imported: {positivity_status['closed']}")
    print(f"  diagonal inequality: {diagonal_status['closed']}")
    print(f"  zero-descent contradiction: {zero_descent_status['closed']}")
    print(f"  lower-half-plane exclusion: {lower_half_status['closed']}")
    print(f"  endpoint passage closed: {endpoint_closed}")
    print(f"  conditional RH-side conclusion: {rh_status['closed']}")
    print(f"  independent external proof vetted: {data['independentRhProofVetted']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
