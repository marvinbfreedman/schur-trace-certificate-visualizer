#!/usr/bin/env python3
r"""Terminal consequence for the Weyl/Volterra external equivalence audit.

The RH/de Branges bridge ledger only needs the external-audit summary:

    rhFacingChainClosed,
    originalKlmConditionClosed,
    normalizedSchurCertificateClosed,
    openCount / blockingItems.

This wrapper reads the full audit and exports only those summary facts, keeping
the detailed evidence list below the full external audit object.
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
        "--external-audit-json",
        default="weyl_volterra_external_equivalence_audit.json",
    )
    parser.add_argument(
        "--json-out",
        default="weyl_volterra_external_audit_consequence_theorem.json",
    )
    args = parser.parse_args()

    audit = load(args.external_audit_json)
    rh_chain = bool(audit.get("rhFacingChainClosed"))
    klm = bool(audit.get("originalKlmConditionClosed"))
    schur = bool(audit.get("normalizedSchurCertificateClosed"))
    external = bool(audit.get("externalFoundationClosed"))
    open_count = int(audit.get("openCount", 0))
    closed_count = int(audit.get("closedCount", audit.get("closedAuditItems", 0) or 0))
    blockers = audit.get("blockingItems", [])
    theorem_closed = rh_chain and klm and schur and external and open_count == 0

    data: dict[str, Any] = {
        "theoremName": "Weyl/Volterra external audit consequence theorem",
        "proofClass": "symbolic identity",
        "source": "Weyl/Volterra external equivalence audit",
        "normalizedSchurCertificateClosed": schur,
        "conditionalQuotientToOriginalLiftClosed": bool(
            audit.get("conditionalQuotientToOriginalLiftClosed")
        ),
        "externalFoundationClosed": external,
        "originalWeylKernelPositivityClosed": bool(
            audit.get("originalWeylKernelPositivityClosed")
        ),
        "originalKlmConditionClosed": klm,
        "rhFacingChainClosed": rh_chain,
        "closedCount": closed_count,
        "openCount": open_count,
        "blockingItems": blockers,
        "nextProofTarget": audit.get("nextProofTarget"),
        "statuses": {
            "normalizedSchurCertificateStatus": status(
                "normalized Schur certificate summary",
                schur,
                "The full audit records the normalized Weyl/Volterra Schur certificate as closed.",
            ),
            "originalKlmConditionStatus": status(
                "original KLM condition summary",
                klm,
                "The full audit records the original KLM condition as closed.",
            ),
            "rhFacingChainStatus": status(
                "RH-facing chain summary",
                rh_chain,
                "The full audit records the RH-facing chain as closed.",
            ),
            "externalAuditConsequenceStatus": status(
                "external audit consequence",
                theorem_closed,
                "Only the audit summary facts are exported upstream.",
            ),
        },
        "externalAuditConsequenceClosed": theorem_closed,
        "proof": [
            "Import the full Weyl/Volterra external equivalence audit.",
            "Export only closed/open counts, blocking items, and terminal chain status flags.",
        ],
        "notExportedHere": [
            "audit item evidence list",
            "quotient theorem filename",
            "bridge theorem filename",
            "source and high-block evidence filenames",
        ],
        "remainingAnalyticGap": None
        if theorem_closed
        else "External audit has open items or missing terminal status flags.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Weyl/Volterra external audit consequence theorem")
    print(f"  normalized Schur: {schur}")
    print(f"  original KLM: {klm}")
    print(f"  RH-facing chain: {rh_chain}")
    print(f"  open count: {open_count}")
    print(f"  theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
