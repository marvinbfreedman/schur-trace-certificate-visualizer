#!/usr/bin/env python3
r"""Compatibility wrapper for the adjoint Green endpoint range theorem.

The old endpoint range theorem routed through the continuum active trace range
criterion.  That was the right fallback when endpoint full row rank was not
certified.  The current proof route is sharper: the synchronized interval /
Krawczyk endpoint theorem proves that the exact active endpoint map has full
active row rank.  Hence the endpoint Fredholm compatibility condition is
vacuous for active source rows.

This wrapper preserves the historical JSON interface while importing the
closed interval theorem:

    adjoint_green_endpoint_range_interval_theorem.json

No sampled-flow endpoint diagnostic is used as proof input here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {"label": label, "closed": ok, "status": "closed" if ok else "open", "reason": reason}


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--interval-endpoint-json",
        default="adjoint_green_endpoint_range_interval_theorem.json",
    )
    # Legacy args are accepted so older scripts do not break, but they are not
    # proof dependencies of this wrapper.
    parser.add_argument("--endpoint-selection-json", default="")
    parser.add_argument("--jump-json", default="")
    parser.add_argument("--active-range-json", default="")
    parser.add_argument("--continuum-active-range-json", default="")
    parser.add_argument("--json-out", default="adjoint_green_endpoint_range_theorem.json")
    args = parser.parse_args()

    interval = load(args.interval_endpoint_json)

    exact_map_ok = closed(interval, "exactFundamentalMatrixEndpointMapStatus")
    fredholm_ok = closed(interval, "endpointFredholmRangeAlternativeStatus")
    full_rank_ok = closed(interval, "endpointFullActiveRowRankStatus")
    bvp_ok = closed(interval, "endpointGreenBvpSolvabilityStatus")
    bound_ok = closed(interval, "uniformTraceDualBoundFromCompactnessStatus")
    compat_ok = closed(interval, "continuumActiveRangeCompatibilityStatus")
    theorem_ok = bool(interval.get("theoremClosed")) and all(
        [exact_map_ok, fredholm_ok, full_rank_ok, bvp_ok, bound_ok, compat_ok]
    )

    data = {
        "theoremName": "exact adjoint Green endpoint range theorem",
        "proofClass": "interval/ball certificate",
        "intervalEndpointJson": args.interval_endpoint_json,
        "legacyInputs": {
            "endpointSelectionJson": args.endpoint_selection_json or None,
            "jumpJson": args.jump_json or None,
            "activeRangeJson": args.active_range_json or None,
            "continuumActiveRangeJson": args.continuum_active_range_json or None,
        },
        "endpointIntervalMargins": interval.get("endpointIntervalMargins", {}),
        "exactFundamentalMatrixEndpointMapStatus": status(
            "exact fundamental-matrix endpoint map",
            exact_map_ok,
            "Imported from the synchronized interval endpoint theorem.",
        ),
        "endpointFredholmRangeAlternativeStatus": status(
            "endpoint Fredholm range alternative",
            fredholm_ok,
            "Imported from the synchronized interval endpoint theorem.",
        ),
        "endpointFullActiveRowRankStatus": status(
            "endpoint map has full active row rank",
            full_rank_ok,
            (
                "The interval/Krawczyk certificate keeps the exact endpoint map "
                "inside the persistent Pluecker chart, so it is onto the active "
                "endpoint row space."
            ),
        ),
        "endpointGreenBvpSolvabilityStatus": status(
            "endpoint Green BVP solvable for every active source row",
            bvp_ok,
            (
                "Because the active endpoint map is onto, Mz=-b(d) is solvable "
                "for every active source row.  No separate continuum active "
                "trace range compatibility is needed at this endpoint layer."
            ),
        ),
        "uniformTraceDualBoundFromCompactnessStatus": status(
            "uniform trace-dual bound on Green solution family",
            bound_ok,
            "Imported from the interval theorem: bounded right inverse plus compact source window.",
        ),
        "continuumActiveRangeCompatibilityStatus": status(
            "continuum active endpoint compatibility",
            compat_ok,
            "Compatibility is automatic because the active endpoint map has full row rank.",
        ),
        "exactEndpointReductionClosed": exact_map_ok and fredholm_ok,
        "unconditionalContinuumEndpointClosed": theorem_ok,
        "theoremClosed": theorem_ok,
        "definitions": interval.get("definitions", {}),
        "formalProof": [
            "Import the synchronized interval endpoint theorem.",
            "Use its exact fundamental-matrix endpoint map and Fredholm alternative.",
            "Use its full active row-rank certificate to make the endpoint compatibility condition vacuous.",
            "Use the bounded right inverse and compact source window for uniform trace-dual bounds.",
        ],
        "remainingAnalyticGap": None if theorem_ok else "The interval endpoint theorem is not closed.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Exact adjoint Green endpoint range theorem")
    print(f"  interval endpoint input: {bool(interval.get('theoremClosed'))}")
    print(f"  exact endpoint reduction: {data['exactEndpointReductionClosed']}")
    print(f"  unconditional endpoint closed: {theorem_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
