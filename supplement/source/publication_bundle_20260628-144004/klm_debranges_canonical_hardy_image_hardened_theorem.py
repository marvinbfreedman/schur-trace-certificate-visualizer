#!/usr/bin/env python3
r"""Hardened canonical Hardy image theorem.

This theorem object separates the exact Hardy/de Branges identity from the
older probe script.  No quadrature or sampled evaluation set is used here.

For an entire Hermite--Biehler candidate E, define

    K_E(w,z) = (E(z)conj(E(w)) - E#(z)conj(E#(w)))
               /(2 pi i (conj(w)-z)).

For Im z>0 and Im w>0, set

    h_z^+(r) = (2 pi)^(-1/2) E(z)  exp(i z r),
    h_z^-(r) = (2 pi)^(-1/2) E#(z) exp(i z r),
    r >= 0.

The identity

    int_0^infty exp(i(z-conj(w))r) dr = 1/(i(conj(w)-z))

is absolutely convergent because Im z+Im w>0.  Therefore

    K_E(w,z) = <h_z^+,h_w^+> - <h_z^-,h_w^->_{L^2(0,infty)}.

For the shifted Riemann Xi family used in the ledger,

    E_omega(z)=Xi(z+i omega),   E_omega#(z)=Xi(z-i omega),

this gives the canonical Hardy branch image of the shifted-Xi de Branges
kernel as a symbolic identity.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json-out",
        default="klm_debranges_canonical_hardy_image_hardened_theorem.json",
    )
    args = parser.parse_args()

    integral_status = status(
        "half-line Hardy integral",
        True,
        (
            "For Im z, Im w > 0, exp(i(z-conj(w))r) decays like "
            "exp(-(Im z+Im w)r).  Direct integration gives "
            "int_0^infty exp(i(z-conj(w))r)dr=1/(i(conj(w)-z))."
        ),
    )
    plus_status = status(
        "plus branch Gram identity",
        True,
        (
            "Substituting h_z^+(r)=(2*pi)^(-1/2)E(z)exp(i z r) into the "
            "L2(0,infty) inner product gives "
            "E(z)conj(E(w))/(2*pi*i*(conj(w)-z))."
        ),
    )
    minus_status = status(
        "minus branch Gram identity",
        True,
        (
            "The same calculation with E# gives "
            "E#(z)conj(E#(w))/(2*pi*i*(conj(w)-z))."
        ),
    )
    kernel_status = status(
        "canonical Hardy/de Branges kernel identity",
        True,
        (
            "Subtracting the minus branch Gram from the plus branch Gram "
            "is exactly the de Branges kernel formula."
        ),
    )

    data = {
        "theoremName": "hardened canonical Hardy image theorem",
        "proofClass": "symbolic identity",
        "domain": {
            "z": "Im z > 0",
            "w": "Im w > 0",
            "r": "0 <= r < infinity",
        },
        "shiftedXiSpecialization": {
            "Eomega": "E_omega(z)=Xi(z+i omega)",
            "EomegaSharp": "E_omega#(z)=Xi(z-i omega)",
            "omegaRangeUsedHere": "real omega; no omega estimate is needed for this identity",
        },
        "formulas": {
            "kernel": (
                "K_E(w,z)=(E(z)conj(E(w))-E#(z)conj(E#(w)))"
                "/(2*pi*i*(conj(w)-z))"
            ),
            "hPlus": "h_z^+(r)=(2*pi)^(-1/2)E(z)exp(i z r)",
            "hMinus": "h_z^-(r)=(2*pi)^(-1/2)E#(z)exp(i z r)",
            "hardyIntegral": (
                "int_0^infty exp(i(z-conj(w))r)dr=1/(i*(conj(w)-z))"
            ),
            "branchIdentity": (
                "K_E(w,z)=<h_z^+,h_w^+>_{L2(0,infty)}"
                "-<h_z^-,h_w^->_{L2(0,infty)}"
            ),
        },
        "statuses": {
            "halfLineHardyIntegralStatus": integral_status,
            "plusBranchGramStatus": plus_status,
            "minusBranchGramStatus": minus_status,
            "canonicalHardyImageStatus": kernel_status,
        },
        "canonicalHardyImageClosed": True,
        "canonicalHardyIdentityClosed": True,
        "symbolicIdentityClosed": True,
        "proof": [
            (
                "Because Im(z-conj(w))=Im z+Im w>0, the exponential "
                "exp(i(z-conj(w))r) is integrable on the positive half-line."
            ),
            (
                "Evaluating the antiderivative gives the Hardy denominator "
                "1/(i(conj(w)-z))."
            ),
            (
                "The plus branch inner product equals the E(z)conj(E(w)) "
                "term of K_E; the minus branch inner product equals the "
                "E#(z)conj(E#(w)) term."
            ),
            "Their signed difference is precisely K_E(w,z).",
        ],
        "scope": (
            "This theorem proves the canonical Hardy branch image only.  The "
            "separate augmented closed-cone theorem supplies the KLM pullback "
            "and positive-cone closure layers."
        ),
        "remainingAnalyticGap": None,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Hardened canonical Hardy image theorem")
    print(f"  half-line integral: {integral_status['closed']}")
    print(f"  plus branch: {plus_status['closed']}")
    print(f"  minus branch: {minus_status['closed']}")
    print(f"  canonical identity: {kernel_status['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
