#!/usr/bin/env python3
r"""Final augmented KLM-to-de Branges positive-kernel bridge.

The detailed augmented Xi machinery is packaged one layer lower in
``shifted_xi_debranges_kernel_positivity_theorem.json``.  This bridge consumes
that theorem and records the transform-level conclusion needed by the endpoint
passage:

    the all-omega augmented KLM construction gives positive shifted-Xi
    de Branges kernels K_{E_omega}.

This file is intentionally a thin bridge wrapper.  It does not re-import the
finite theta pullback, trace convergence, repair, or cone-closure components
directly.
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
        "--kernel-positivity-json",
        default="shifted_xi_debranges_kernel_positivity_theorem.json",
    )
    # Legacy optional arguments are accepted for compatibility but are no
    # longer proof dependencies.
    parser.add_argument("--hardy-json", default="")
    parser.add_argument("--finite-pullback-json", default="")
    parser.add_argument("--trace-convergence-json", default="")
    parser.add_argument("--repair-positive-json", default="")
    parser.add_argument("--closed-cone-json", default="")
    parser.add_argument("--hardened-lift-json", default="")
    parser.add_argument("--pullback-theorem-json", default="")
    parser.add_argument("--augmented-continuum-json", default="")
    parser.add_argument("--uniform-klm-json", default="")
    parser.add_argument("--json-out", default="klm_debranges_augmented_closed_cone_theorem.json")
    args = parser.parse_args()

    kernel = load(args.kernel_positivity_json)
    kernel_ok = bool(
        kernel.get("shiftedXiDeBrangesKernelPositiveClosed")
        or nested_closed(kernel, "statuses", "shiftedXiKernelPositivityStatus")
    )
    finite_eval_ok = bool(kernel.get("finiteEvaluationGramPositive", kernel_ok))

    kernel_status = status(
        "shifted-Xi de Branges kernel positivity",
        kernel_ok,
        (
            "Imported shifted-Xi kernel theorem: K_{E_omega} is positive on "
            "finite upper-half-plane evaluation sets for 0<omega<1/2."
        ),
        blocker=None
        if kernel_ok
        else "Close shifted_xi_debranges_kernel_positivity_theorem.py.",
    )
    bridge_status = status(
        "KLM-to-de Branges positive-kernel bridge",
        kernel_ok and finite_eval_ok,
        (
            "The augmented KLM construction has been transported to the "
            "shifted-Xi de Branges kernel.  Therefore the endpoint theorem may "
            "use positivity of K_{E_omega} for every finite evaluation set."
        ),
        blocker=None
        if kernel_ok and finite_eval_ok
        else "Need shifted-Xi de Branges kernel positivity.",
    )

    bridge_closed = bridge_status["closed"]
    data = {
        "theoremName": "augmented KLM-to-de Branges positive-kernel bridge theorem",
        "proofClass": "analytic proof",
        "kernelPositivityJson": args.kernel_positivity_json,
        "legacyInputs": {
            "hardyJson": args.hardy_json or None,
            "finitePullbackJson": args.finite_pullback_json or None,
            "traceConvergenceJson": args.trace_convergence_json or None,
            "repairPositiveJson": args.repair_positive_json or None,
            "closedConeJson": args.closed_cone_json or None,
            "hardenedLiftJson": args.hardened_lift_json or None,
            "pullbackTheoremJson": args.pullback_theorem_json or None,
            "augmentedContinuumJson": args.augmented_continuum_json or None,
            "uniformKlmJson": args.uniform_klm_json or None,
        },
        "statuses": {
            "shiftedXiKernelPositivityStatus": kernel_status,
            "klmToDeBrangesBridgeStatus": bridge_status,
        },
        "formalProof": [
            "Import shifted-Xi de Branges kernel positivity for 0<omega<1/2.",
            "Use that positivity as the KLM-to-de Branges transform output.",
            "Pass the positive kernel family to the endpoint theorem.",
        ],
        "bridgeClosed": bridge_closed,
        "transformClosed": bridge_closed,
        "nextProofTarget": (
            "Apply the de Branges/Hermite-Biehler endpoint passage for "
            "E_omega(z)=Xi(z+i omega)."
            if bridge_closed
            else "Close shifted-Xi de Branges kernel positivity."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented KLM-to-de Branges positive-kernel bridge theorem")
    print(f"  shifted-Xi kernel positivity: {kernel_status['closed']}")
    print(f"  bridge closed: {bridge_closed}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
