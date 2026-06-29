#!/usr/bin/env python3
r"""Zero-descent endpoint theorem for shifted Xi."""

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
        "--diagonal-json",
        default="debranges_diagonal_inequality_theorem.json",
    )
    parser.add_argument(
        "--normalization-json",
        default="shifted_xi_debranges_normalization_theorem.json",
    )
    parser.add_argument("--json-out", default="shifted_xi_zero_descent_endpoint_theorem.json")
    args = parser.parse_args()

    diagonal = load(args.diagonal_json)
    normalization = load(args.normalization_json)

    diagonal_ok = bool(diagonal.get("diagonalInequalityClosed"))
    real_entire_ok = bool(normalization.get("EomegaSharpIdentityClosed"))
    zero_descent_ok = diagonal_ok
    lower_half_ok = zero_descent_ok and real_entire_ok
    endpoint_ok = zero_descent_ok and lower_half_ok

    data = {
        "theoremName": "shifted Xi zero-descent endpoint theorem",
        "proofClass": "analytic proof",
        "diagonalJson": args.diagonal_json,
        "normalizationJson": args.normalization_json,
        "statement": (
            "The diagonal shifted-Xi inequality for every 0<omega<1/2 excludes "
            "all non-real zeros of Xi(z)=xi(1/2+i z)."
        ),
        "statuses": {
            "diagonalInequalityInputStatus": status(
                "diagonal inequality input",
                diagonal_ok,
                "|Xi(z-i omega)| <= |Xi(z+i omega)| on C_+.",
            ),
            "zeroDescentContradictionStatus": status(
                "zero-descent contradiction",
                zero_descent_ok,
                (
                    "If Xi(rho)=0 with Im rho=b>0, then z=rho-i omega lies in C_+ "
                    "for 0<omega<min(b,1/2), and the diagonal inequality forces "
                    "Xi(rho-2i omega)=0 on a continuum, impossible for nonzero entire Xi."
                ),
            ),
            "realEntireConjugationStatus": status(
                "real-entire conjugation exclusion",
                lower_half_ok,
                (
                    "Xi is real entire, so lower-half-plane zeros conjugate to "
                    "upper-half-plane zeros, already excluded."
                ),
            ),
            "endpointZeroLocationStatus": status(
                "endpoint zero-location theorem",
                endpoint_ok,
                "All zeros of Xi in the z-normalization are real.",
            ),
        },
        "zeroDescentEndpointClosed": endpoint_ok,
        "endpointZeroLocationClosed": endpoint_ok,
        "conditionalRhClosed": endpoint_ok,
        "proof": [
            "Assume Xi(rho)=0 with Im rho=b>0.",
            "For 0<omega<min(b,1/2), set z=rho-i omega so z lies in C_+.",
            "The diagonal inequality gives |Xi(rho-2i omega)|<=|Xi(rho)|=0.",
            "The resulting vertical continuum of zeros has accumulation points in C.",
            "A nonzero entire function cannot vanish on such a set; conjugation handles the lower half-plane.",
        ],
        "remainingAnalyticGap": None if endpoint_ok else "Diagonal inequality or real-entire normalization is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Shifted-Xi zero-descent endpoint theorem")
    print(f"  diagonal inequality: {diagonal_ok}")
    print(f"  lower-half exclusion: {lower_half_ok}")
    print(f"  endpoint zero location: {endpoint_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
