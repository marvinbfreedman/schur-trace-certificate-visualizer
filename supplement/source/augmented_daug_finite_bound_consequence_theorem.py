#!/usr/bin/env python3
r"""Narrow finite-bound consequence for the augmented D_aug representation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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
        "--finite-constants-json",
        default="xi_augmented_finite_schur_interval_constants_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="augmented_daug_finite_bound_consequence_theorem.json",
    )
    args = parser.parse_args()

    constants = load(args.finite_constants_json)
    nonnegative_ok = bool(constants.get("finiteRepairOperatorNonnegativeConstantsClosed"))
    bound_ok = bool(constants.get("finiteRepairOperatorBoundConstantsClosed"))
    interval = constants.get("intervalCertificate", {})
    norm_upper = interval.get("dMaxUpper")
    ok = nonnegative_ok and bound_ok and norm_upper is not None

    data = {
        "theoremName": "augmented D_aug finite bound consequence theorem",
        "proofClass": "symbolic identity",
        "finiteConstantsJson": args.finite_constants_json,
        "operatorBound": {
            "normUpper": norm_upper,
            "source": "interval dMaxUpper from the finite Schur repair constants",
        },
        "statuses": {
            "finiteRepairNonnegativeInputStatus": status(
                "finite D_aug,N nonnegative constants",
                nonnegative_ok,
                "The interval constants theorem proves finite D_aug,N is nonnegative up to the certified ball.",
            ),
            "finiteRepairBoundInputStatus": status(
                "finite D_aug,N upper-bound constants",
                bound_ok,
                "The interval constants theorem exports dMaxUpper for the finite repair operator.",
            ),
            "finiteBoundConsequenceStatus": status(
                "finite D_aug bound consequence",
                ok,
                "Only nonnegativity and dMaxUpper are exported to the continuum D_aug representation theorem.",
            ),
        },
        "DaugFiniteBoundConsequenceClosed": ok,
        "finiteRepairOperatorNonnegativeConstantsClosed": nonnegative_ok,
        "finiteRepairOperatorBoundConstantsClosed": bound_ok,
        "DaugOperatorNormUpper": norm_upper,
        "proof": [
            "Import the finite Schur interval constants theorem.",
            "Project away all finite Schur details except finite repair nonnegativity and dMaxUpper.",
        ],
        "notExportedHere": [
            "range residual",
            "Schur complement gap",
            "Mu annihilation",
            "finite repaired form positivity",
        ],
        "remainingAnalyticGap": None if ok else "Finite repair nonnegativity or bound constants are open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented D_aug finite bound consequence theorem")
    print(f"  finite bound consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
