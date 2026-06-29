#!/usr/bin/env python3
r"""Terminal consequence of the 200-point endpoint eigenrow propagation certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {"label": label, "closed": closed, "status": "closed" if closed else "open", "reason": reason}
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--eigenrow-json",
        default="endpoint_eigenrow_interval_propagation_200.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_eigenrow_interval_propagation_200_consequence_theorem.json",
    )
    args = parser.parse_args()

    eigenrow = load(args.eigenrow_json)
    center = bool(eigenrow.get("centerSynchronizationStatus", {}).get("closed"))
    target = bool(eigenrow.get("eigenrowTraceTargetStatus", {}).get("closed"))
    krawczyk = bool(eigenrow.get("eigenrowTaylorKrawczykStatus", {}).get("closed"))
    closed = center and target and krawczyk

    data = {
        "theoremName": "endpoint eigenrow interval propagation 200 consequence theorem",
        "proofClass": "symbolic identity",
        "centerSynchronizationStatus": status(
            "center synchronization",
            center,
            "The 200-point eigenrow propagation certificate uses the synchronized quadrature center.",
        ),
        "eigenrowTraceTargetStatus": status(
            "eigenrow trace target",
            target,
            "The propagated eigenrow derivatives meet the trace-radius target.",
        ),
        "eigenrowTaylorKrawczykStatus": status(
            "eigenrow Taylor/Krawczyk propagation",
            krawczyk,
            "The eigenrow Taylor recurrence is enclosed by the Krawczyk certificate.",
        ),
        "conditionMatrixOrder": eigenrow.get("conditionMatrixOrder"),
        "controlledQuadratureOrder": eigenrow.get("controlledQuadratureOrder"),
        "maxDerivativeEntryRadius": eigenrow.get("maxDerivativeEntryRadius"),
        "maxDerivativeEntryRadiusText": eigenrow.get("maxDerivativeEntryRadiusText"),
        "endpointEigenrowPropagation200Closed": closed,
        "remainingAnalyticGap": None if closed else "Close endpoint_eigenrow_interval_propagation_200.py.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint eigenrow interval propagation 200 consequence theorem")
    print(f"  eigenrow propagation closed: {closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
