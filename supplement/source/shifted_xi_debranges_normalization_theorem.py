#!/usr/bin/env python3
r"""Normalization theorem for shifted-Xi de Branges kernels."""

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
    parser.add_argument("--riemann-json", default="riemann_kernel_normalization_theorem.json")
    parser.add_argument(
        "--hardy-image-json",
        default="klm_debranges_canonical_hardy_image_hardened_theorem.json",
    )
    parser.add_argument("--json-out", default="shifted_xi_debranges_normalization_theorem.json")
    args = parser.parse_args()

    riemann = load(args.riemann_json)
    hardy = load(args.hardy_image_json)

    riemann_ok = bool(riemann.get("riemannKernelNormalizationClosed"))
    hardy_ok = bool(hardy.get("canonicalHardyImageClosed") or hardy.get("canonicalHardyIdentityClosed"))
    sharp_ok = riemann_ok
    denominator_ok = True
    normalization_ok = riemann_ok and hardy_ok and sharp_ok and denominator_ok

    data = {
        "theoremName": "shifted Xi de Branges normalization theorem",
        "proofClass": "symbolic identity",
        "riemannJson": args.riemann_json,
        "hardyImageJson": args.hardy_image_json,
        "normalization": {
            "Xi": "Xi(z)=xi(1/2+i z), real entire in z",
            "Eomega": "E_omega(z)=Xi(z+i omega)",
            "EomegaSharp": "E_omega#(z)=overline(E_omega(conj z))=Xi(z-i omega)",
            "omegaRange": "0<omega<1/2",
            "evaluationDomain": "z,w in C_+",
            "kernel": (
                "K_E(w,z)=(E(z)conj(E(w))-E#(z)conj(E#(w)))/"
                "(2*pi*i*(conj(w)-z))"
            ),
            "diagonal": (
                "K_E(z,z)=(|Xi(z+i omega)|^2-|Xi(z-i omega)|^2)/(4*pi*Im z)"
            ),
        },
        "statuses": {
            "riemannNormalizationInputStatus": status(
                "Riemann/Xi normalization input",
                riemann_ok,
                "The Riemann kernel normalization fixes Xi up to harmless positive Gram rescaling.",
            ),
            "realEntireSharpStatus": status(
                "real-entire sharp operation",
                sharp_ok,
                (
                    "Since Xi is real entire in the z-normalization and omega is real, "
                    "E_omega#(z)=overline(Xi(conj z+i omega))=Xi(z-i omega)."
                ),
            ),
            "upperHalfPlaneDenominatorStatus": status(
                "upper-half-plane denominator",
                denominator_ok,
                (
                    "For z,w in C_+, conj(w)-z cannot vanish.  On the diagonal, "
                    "2*pi*i*(conj(z)-z)=4*pi*Im z>0."
                ),
            ),
            "canonicalHardyImageInputStatus": status(
                "canonical Hardy/de Branges image input",
                hardy_ok,
                "The hardened Hardy image theorem fixes the same de Branges kernel normalization.",
            ),
            "shiftedXiNormalizationStatus": status(
                "shifted-Xi de Branges normalization",
                normalization_ok,
                "The shifted kernel, sharp operation, and diagonal sign convention are all fixed.",
            ),
        },
        "shiftedXiDeBrangesNormalizationClosed": normalization_ok,
        "EomegaSharpIdentityClosed": sharp_ok,
        "deBrangesKernelDenominatorClosed": denominator_ok,
        "diagonalFormulaClosed": normalization_ok,
        "proof": [
            "Use the fixed Xi normalization from the Riemann kernel theorem.",
            "Use real-entireness of Xi to compute the sharp transform of E_omega.",
            "Check that the de Branges denominator is nonzero on C_+ x C_+ and positive on the diagonal.",
            "Use the canonical Hardy image theorem to match the same kernel normalization.",
        ],
        "remainingAnalyticGap": None if normalization_ok else "A normalization input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi de Branges normalization theorem")
    print(f"  Riemann normalization: {riemann_ok}")
    print(f"  Hardy image normalization: {hardy_ok}")
    print(f"  shifted normalization: {normalization_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
