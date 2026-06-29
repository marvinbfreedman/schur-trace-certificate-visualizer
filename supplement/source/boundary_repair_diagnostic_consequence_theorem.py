#!/usr/bin/env python3
r"""Narrow diagnostic consequence of the boundary-repair identity ledger."""

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


def nested_closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--boundary-repair-json",
        default="boundary_repair_identity.json",
    )
    parser.add_argument(
        "--json-out",
        default="boundary_repair_diagnostic_consequence_theorem.json",
    )
    args = parser.parse_args()

    repair = load(args.boundary_repair_json)
    primitive_false = bool(repair.get("primitiveVanishingFalse"))
    alternative_a_closed = nested_closed(repair, "alternativeAResolutionStatus")
    nonunique = nested_closed(repair, "abstractRepairNonuniquenessStatus")
    alternative_b_closed = nested_closed(repair, "alternativeBAsStatedStatus")
    ok = primitive_false and alternative_a_closed and nonunique and alternative_b_closed

    primitive_status = status(
        "primitive vanishing route is false",
        primitive_false and alternative_a_closed,
        "Imported from the boundary-repair identity diagnostic ledger.",
    )
    alternative_a_status = status(
        "alternative A cannot be proved",
        alternative_a_closed,
        "Imported from the boundary-repair identity diagnostic ledger.",
    )
    nonunique_status = status(
        "abstract quotient repair is nonunique",
        nonunique and alternative_b_closed,
        "Imported from the boundary-repair identity diagnostic ledger.",
    )
    alternative_b_status = status(
        "alternative B as stated is not well-defined",
        alternative_b_closed,
        "Imported from the boundary-repair identity diagnostic ledger.",
    )

    data = {
        "theoremName": "boundary repair diagnostic consequence theorem",
        "proofClass": "symbolic identity",
        "boundaryRepairJson": args.boundary_repair_json,
        "statement": (
            "The boundary-repair identity ledger supplies only diagnostic "
            "facts used by the canonical boundary comparison: primitive "
            "vanishing is false and the abstract quotient repair is nonunique. "
            "It does not claim the full boundary-repair identity is closed."
        ),
        "statuses": {
            "primitiveVanishingFalseStatus": primitive_status,
            "alternativeAResolutionStatus": alternative_a_status,
            "abstractRepairNonuniquenessStatus": nonunique_status,
            "alternativeBAsStatedStatus": alternative_b_status,
            "boundaryRepairDiagnosticConsequenceStatus": status(
                "boundary repair diagnostic consequence",
                ok,
                (
                    "The comparison theorem may use these diagnostic facts "
                    "while the full boundary-repair identity remains a separate "
                    "open/retired target."
                ),
            ),
        },
        "primitiveVanishingFalse": primitive_false,
        "primitiveVanishingRouteStatus": primitive_status,
        "alternativeAResolutionStatus": alternative_a_status,
        "abstractRepairNonuniquenessStatus": nonunique_status,
        "alternativeBAsStatedStatus": alternative_b_status,
        "alternativeAFalseClosed": alternative_a_closed,
        "abstractRepairNonuniquenessClosed": nonunique,
        "alternativeBAsStatedNotWellDefinedClosed": alternative_b_closed,
        "boundaryRepairDiagnosticConsequenceClosed": ok,
        "boundaryRepairIdentityClosed": bool(repair.get("boundaryRepairIdentityClosed")),
        "proof": [
            "Import the boundary-repair identity diagnostic ledger.",
            "Extract the primitive-vanishing disproof.",
            "Extract the abstract quotient-repair nonuniqueness statement.",
            "Do not import the intentionally open full boundary-repair identity as a closed proof.",
        ],
        "remainingAnalyticGap": None
        if ok
        else "Boundary diagnostic facts are not all closed.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Boundary repair diagnostic consequence theorem")
    print(f"  primitive vanishing false: {primitive_false and alternative_a_closed}")
    print(f"  abstract repair nonunique: {nonunique and alternative_b_closed}")
    print(f"  consequence: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
