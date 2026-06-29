#!/usr/bin/env python3
r"""Terminal consequence for the uniform-omega Weyl/KLM bridge.

The RH/de Branges bridge ledger only needs the all-omega KLM/Weyl positivity
facts:

    uniformOmegaCoverageClosed = true,
    originalKlmConditionClosed = true.

This wrapper reads the full uniform-omega Weyl/KLM bridge and exports only
those terminal facts, keeping the Weyl operator-family and hbar=1 KLM/Weyl
equivalence details below the full bridge theorem.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: str) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict[str, Any]:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--uniform-json",
        default="uniform_omega_weyl_klm_bridge.json",
    )
    parser.add_argument(
        "--json-out",
        default="uniform_omega_weyl_klm_consequence_theorem.json",
    )
    args = parser.parse_args()

    uniform = load(args.uniform_json)
    omega_ok = bool(uniform.get("uniformOmegaCoverageClosed"))
    weyl_ok = bool(uniform.get("originalWeylOperatorPositiveClosed"))
    klm_ok = bool(uniform.get("originalKlmConditionClosed"))
    theorem_closed = omega_ok and weyl_ok and klm_ok

    data: dict[str, Any] = {
        "theoremName": "uniform omega Weyl/KLM consequence theorem",
        "proofClass": "symbolic identity",
        "source": "uniform omega Weyl/KLM bridge",
        "omegaRange": uniform.get("omegaRange", "|omega| < 1/2"),
        "statuses": {
            "uniformOmegaCoverageInputStatus": status(
                "uniform omega coverage input",
                omega_ok,
                "The full bridge proves coverage for every |omega|<1/2.",
            ),
            "weylOperatorPositivityInputStatus": status(
                "positive Weyl operator family input",
                weyl_ok,
                "The full bridge imports the positive Weyl operator family.",
            ),
            "klmPositiveTypeInputStatus": status(
                "KLM positive-type input",
                klm_ok,
                "The full bridge applies the hbar=1 KLM/Weyl equivalence.",
            ),
            "uniformOmegaKlmConsequenceStatus": status(
                "uniform omega KLM/Weyl consequence",
                theorem_closed,
                "Only the terminal all-omega KLM/Weyl positivity facts are exported upstream.",
            ),
        },
        "uniformOmegaCoverageClosed": omega_ok,
        "originalWeylKernelPositivityClosed": bool(
            uniform.get("originalWeylKernelPositivityClosed")
        ),
        "originalWeylOperatorPositiveClosed": weyl_ok,
        "originalKlmConditionClosed": klm_ok,
        "uniformOmegaWeylKlmConsequenceClosed": theorem_closed,
        "proof": [
            "Import the full uniform-omega Weyl/KLM bridge.",
            "Export only all-omega coverage and KLM positive-type closure for the RH/de Branges ledger.",
        ],
        "notExportedHere": [
            "positive Weyl operator-family proof",
            "hbar=1 KLM/Weyl equivalence proof",
            "normalization details",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "Uniform-omega Weyl/KLM bridge is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Uniform omega Weyl/KLM consequence theorem")
    print(f"  uniform omega coverage: {omega_ok}")
    print(f"  Weyl operator positivity: {weyl_ok}")
    print(f"  KLM positive type: {klm_ok}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
