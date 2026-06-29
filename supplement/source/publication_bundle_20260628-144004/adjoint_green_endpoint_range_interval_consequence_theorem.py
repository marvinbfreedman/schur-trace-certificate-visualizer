#!/usr/bin/env python3
r"""Terminal consequence of the adjoint Green endpoint range interval theorem."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def item_closed(data: dict, key: str) -> bool:
    return bool(data.get(key, {}).get("closed"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {"label": label, "closed": closed, "status": "closed" if closed else "open", "reason": reason}
    if blocker:
        out["blocker"] = blocker
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--endpoint-range-json",
        default="adjoint_green_endpoint_range_interval_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="adjoint_green_endpoint_range_interval_consequence_theorem.json",
    )
    args = parser.parse_args()

    endpoint = load(args.endpoint_range_json)
    exact = item_closed(endpoint, "exactFundamentalMatrixEndpointMapStatus")
    fredholm = item_closed(endpoint, "endpointFredholmRangeAlternativeStatus")
    bvp = item_closed(endpoint, "endpointGreenBvpSolvabilityStatus")
    bound = item_closed(endpoint, "uniformTraceDualBoundFromCompactnessStatus")
    theorem_closed = bool(endpoint.get("theoremClosed")) and exact and fredholm and bvp and bound

    data = {
        "theoremName": "adjoint Green endpoint range interval consequence theorem",
        "proofClass": "symbolic identity",
        "exactFundamentalMatrixEndpointMapStatus": status(
            "exact fundamental-matrix endpoint map",
            exact,
            "The endpoint range theorem reduces active endpoint compatibility to the finite endpoint map.",
        ),
        "endpointFredholmRangeAlternativeStatus": status(
            "endpoint Fredholm range alternative",
            fredholm,
            "The endpoint range theorem identifies solvability with the finite Fredholm range condition.",
        ),
        "endpointGreenBvpSolvabilityStatus": status(
            "endpoint Green BVP solvability",
            bvp,
            "The synchronized endpoint interval certificate makes the active endpoint Fredholm obstruction vacuous.",
        ),
        "uniformTraceDualBoundFromCompactnessStatus": status(
            "uniform trace-dual endpoint bound",
            bound,
            "The endpoint range theorem supplies the uniform trace-dual bound needed for closure.",
        ),
        "theoremClosed": theorem_closed,
        "unconditionalContinuumEndpointClosed": bool(endpoint.get("unconditionalContinuumEndpointClosed")),
        "remainingAnalyticGap": None if theorem_closed else "Close adjoint_green_endpoint_range_interval_theorem.py.",
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Adjoint Green endpoint range interval consequence theorem")
    print(f"  endpoint range theorem closed: {theorem_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
