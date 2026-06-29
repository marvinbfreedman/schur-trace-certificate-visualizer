#!/usr/bin/env python3
r"""Narrow quotient minimal-repair consequence theorem.

This exports only the quotient-Schur facts needed by the canonical boundary
repair comparison:

    Q_Phi = ||G_q f||^2 - <D_q Rf,Rf>
    D_q = (Gamma^*Gamma-C)_+ is bounded on X_R.

The active/inactive source proof and high-block exhaustion details remain in
``weyl_volterra_quotient_schur_theorem.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else bool(item)


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
        "--quotient-json",
        default="weyl_volterra_quotient_schur_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="quotient_minimal_repair_consequence_theorem.json",
    )
    args = parser.parse_args()

    quotient = load(args.quotient_json)
    schur_ok = closed(quotient, "globalWeylVolterraSchurStatus")
    factorization_ok = closed(quotient, "quotientFactorizationStatus")
    repair_ok = closed(quotient, "traceSideRepairStatus")
    moore_ok = closed(quotient, "moorePenroseSchurStatus")
    ok = schur_ok and factorization_ok and repair_ok and moore_ok

    data = {
        "theoremName": "quotient minimal repair consequence theorem",
        "proofClass": "symbolic identity",
        "transportedTraceSpace": {
            "name": "X_R",
            "definition": "completed transported trace range R(U) with quotient norm",
        },
        "statuses": {
            "quotientSchurInputStatus": status(
                "normalized quotient Schur input",
                schur_ok,
                "The quotient Schur theorem closes the normalized full-Phi quotient model.",
            ),
            "quotientFactorizationInputStatus": status(
                "quotient factorization input",
                factorization_ok,
                "The quotient theorem gives Q_Phi=||G_q f||^2-||S Rf||^2.",
            ),
            "traceSideRepairInputStatus": status(
                "bounded trace-side repair input",
                repair_ok,
                "The quotient theorem exports a bounded trace-side repair operator on X_R.",
            ),
            "minimalRepairInputStatus": status(
                "Moore-Penrose minimal repair input",
                moore_ok,
                "The quotient theorem identifies the singular Moore-Penrose Schur repair.",
            ),
            "minimalRepairConsequenceStatus": status(
                "minimal repair consequence",
                ok,
                "The imported quotient facts define D_q=(Gamma^*Gamma-C)_+ as a bounded repair on X_R.",
            ),
        },
        "quotientMinimalRepairConsequenceClosed": ok,
        "minimalQuotientRepairClosed": ok,
        "quotientFactorizationClosed": factorization_ok,
        "traceSideRepairClosed": repair_ok,
        "globalWeylVolterraSchurClosed": schur_ok,
        "operator": {
            "name": "D_q",
            "definition": "(Gamma^*Gamma-C)_+",
            "space": "X_R",
            "claim": "bounded minimal trace repair for the fixed quotient positive part",
        },
        "proof": [
            "Import the quotient Schur theorem.",
            "Expose only the factorization and bounded minimal-repair consequences.",
        ],
        "notExportedHere": [
            "active range proof",
            "source-inactive high-block proof",
            "finite/source quadrature details",
        ],
        "remainingAnalyticGap": None if ok else "One quotient Schur minimal-repair input is open.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Quotient minimal repair consequence theorem")
    print(f"  minimal repair consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
