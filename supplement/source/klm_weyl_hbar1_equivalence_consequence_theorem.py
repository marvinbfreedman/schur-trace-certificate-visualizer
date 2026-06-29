#!/usr/bin/env python3
r"""Terminal consequence of the hbar=1 KLM/Weyl equivalence theorem.

The uniform-omega bridge does not need the full convention theorem.  It only
needs the closed consequence

    KLM positive type of Q <=> positivity of Op^W(sigma)

in the fixed hbar=1 symplectic Fourier normalization.  This wrapper exports
only that terminal fact so upstream proof ledgers do not depend directly on the
full normalization audit node.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--equivalence-json",
        default="klm_weyl_hbar1_equivalence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="klm_weyl_hbar1_equivalence_consequence_theorem.json",
    )
    args = parser.parse_args()

    theorem = load(args.equivalence_json)
    closed = bool(theorem.get("klmWeylHbar1EquivalenceClosed"))

    consequence_status = status(
        "hbar=1 KLM/Weyl equivalence consequence",
        closed,
        (
            "The fixed hbar=1 symplectic Fourier convention identifies KLM "
            "quantum positive type of Q with positivity of the Weyl operator "
            "Op^W(sigma), for the paired objects used in the ledger."
        ),
        blocker=None if closed else "Close klm_weyl_hbar1_equivalence_theorem.py.",
    )

    data = {
        "theoremName": "hbar=1 KLM/Weyl equivalence consequence theorem",
        "proofClass": "symbolic identity",
        "normalization": theorem.get("normalization", "hbar=1"),
        "statuses": {
            "klmWeylEquivalenceConsequenceStatus": consequence_status,
        },
        "klmWeylHbar1EquivalenceClosed": consequence_status["closed"],
        "klmPositiveTypeEqualsWeylOperatorPositiveClosed": consequence_status["closed"],
        "proof": [
            "Import the full hbar=1 KLM/Weyl convention theorem.",
            "Export only the terminal equivalence consequence needed by the uniform-omega bridge.",
        ],
        "remainingAnalyticGap": None if consequence_status["closed"] else consequence_status["blocker"],
        "nextProofTarget": None if consequence_status["closed"] else consequence_status["blocker"],
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("hbar=1 KLM/Weyl equivalence consequence theorem")
    print(f"  consequence closed: {data['klmWeylHbar1EquivalenceClosed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
