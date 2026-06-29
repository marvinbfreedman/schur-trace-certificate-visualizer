#!/usr/bin/env python3
r"""Diagonal de Branges inequality for shifted Xi."""

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
        "--kernel-positivity-json",
        default="shifted_xi_debranges_kernel_positivity_theorem.json",
    )
    parser.add_argument(
        "--normalization-json",
        default="shifted_xi_debranges_normalization_theorem.json",
    )
    parser.add_argument("--json-out", default="debranges_diagonal_inequality_theorem.json")
    args = parser.parse_args()

    kernel = load(args.kernel_positivity_json)
    normalization = load(args.normalization_json)

    kernel_ok = bool(kernel.get("shiftedXiDeBrangesKernelPositiveClosed"))
    diagonal_ok = bool(normalization.get("diagonalFormulaClosed"))
    inequality_ok = kernel_ok and diagonal_ok

    data = {
        "theoremName": "shifted Xi diagonal de Branges inequality theorem",
        "proofClass": "symbolic identity",
        "kernelPositivityJson": args.kernel_positivity_json,
        "normalizationJson": args.normalization_json,
        "statement": (
            "For z in C_+ and 0<omega<1/2, shifted-Xi de Branges kernel "
            "positivity implies |Xi(z-i omega)| <= |Xi(z+i omega)|."
        ),
        "statuses": {
            "kernelPositivityInputStatus": status(
                "shifted kernel positivity input",
                kernel_ok,
                "K_{E_omega}(z,z)>=0 for every z in C_+.",
            ),
            "diagonalFormulaInputStatus": status(
                "diagonal formula input",
                diagonal_ok,
                "K_E(z,z) has denominator 4*pi*Im z and the shifted-Xi numerator.",
            ),
            "diagonalInequalityStatus": status(
                "diagonal shifted-Xi inequality",
                inequality_ok,
                "Since 4*pi*Im z>0, K_E(z,z)>=0 gives the stated modulus inequality.",
            ),
        },
        "diagonalInequalityClosed": inequality_ok,
        "shiftedXiDiagonalInequalityClosed": inequality_ok,
        "proof": [
            "Apply positive definiteness to the one-point evaluation set {z}.",
            "Use the diagonal normalization formula.",
            "Multiply by the positive denominator 4*pi*Im z.",
        ],
        "remainingAnalyticGap": None if inequality_ok else "Kernel positivity or diagonal normalization is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi diagonal inequality theorem")
    print(f"  kernel positivity: {kernel_ok}")
    print(f"  diagonal formula: {diagonal_ok}")
    print(f"  diagonal inequality: {inequality_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
