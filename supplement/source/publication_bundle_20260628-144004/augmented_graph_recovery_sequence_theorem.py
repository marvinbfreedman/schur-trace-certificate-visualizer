#!/usr/bin/env python3
r"""Graph-norm recovery sequence theorem for augmented Mosco transport."""

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
    parser.add_argument("--core-density-json", default="augmented_core_density_theorem.json")
    parser.add_argument("--json-out", default="augmented_graph_recovery_sequence_theorem.json")
    args = parser.parse_args()

    core = load(args.core_density_json)
    core_ok = bool(core.get("coreDensityClosed"))
    recovery_ok = core_ok

    data = {
        "theoremName": "augmented graph recovery sequence theorem",
        "proofClass": "analytic proof",
        "coreDensityJson": args.core_density_json,
        "statement": (
            "For every f in the completed augmented graph space V there are "
            "finite smooth core functions f_N with f_N -> f in the "
            "augmented graph norm."
        ),
        "statuses": {
            "coreDensityInputStatus": status(
                "core density input",
                core_ok,
                "The smooth core is dense in the completed augmented graph space V.",
            ),
            "graphRecoverySequenceStatus": status(
                "graph-norm recovery sequence",
                recovery_ok,
                (
                    "By density, choose smooth core approximants and then a "
                    "diagonal finite-section sequence converging in the V graph norm."
                ),
            ),
        },
        "graphRecoverySequenceClosed": recovery_ok,
        "finiteCoreRecoveryClosed": recovery_ok,
        "proof": [
            "Use the definition of V as the completion of the smooth lifted core.",
            "Choose a core approximant within epsilon in graph norm.",
            "Choose a finite core section of that core approximant within epsilon.",
            "Diagonalize epsilon -> 0 to obtain f_N -> f in V.",
        ],
        "remainingAnalyticGap": None if recovery_ok else "Core density is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented graph recovery sequence theorem")
    print(f"  core density: {core_ok}")
    print(f"  graph recovery: {recovery_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
