#!/usr/bin/env python3
r"""Narrow endpoint full-active-row-rank consequence theorem.

This wrapper exposes only the consequence needed by the synchronized active
range theorem:

    exact endpoint map has full active row rank
    => endpoint map is onto the active row space
    => a bounded endpoint right inverse exists.

The detailed interval/Krawczyk coefficient and boundary-row bookkeeping stays
inside ``endpoint_coefficient_synchronized_200_certificate.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--endpoint-json",
        default="endpoint_coefficient_synchronized_200_certificate.json",
    )
    parser.add_argument(
        "--json-out",
        default="endpoint_full_active_row_rank_consequence_theorem.json",
    )
    args = parser.parse_args()

    endpoint = load(args.endpoint_json)

    center_closed = closed(endpoint, "centerSynchronizationStatus")
    companion_closed = closed(endpoint, "analyticCompanionInputStatus")
    boundary_closed = closed(endpoint, "analyticBoundaryInputStatus")
    krawczyk_closed = closed(endpoint, "synchronizedKrawczykInputStatus")

    q = float(endpoint.get("actualKrawczykQ", float("inf")))
    companion_ratio = float(
        endpoint.get("scaledCompanionRadiusOverSimultaneousBudget", float("inf"))
    )
    boundary_ratio = float(
        endpoint.get("boundaryRowRadiusOverSimultaneousBudget", float("inf"))
    )
    chart_ratio = float(endpoint.get("endpointRadiusOverChartCapacity", float("inf")))

    margins_ok = (
        q < 1.0
        and companion_ratio < 1.0
        and boundary_ratio < 1.0
        and chart_ratio < 1.0
    )
    endpoint_input_ok = (
        center_closed
        and companion_closed
        and boundary_closed
        and krawczyk_closed
        and margins_ok
    )
    pluecker_ok = endpoint_input_ok
    full_rank_ok = pluecker_ok
    onto_ok = full_rank_ok
    right_inverse_ok = onto_ok

    data = {
        "theoremName": "endpoint full active row rank consequence theorem",
        "proofClass": "symbolic identity",
        "conditionMatrixOrder": endpoint.get("conditionMatrixOrder"),
        "controlledQuadratureOrder": endpoint.get("controlledQuadratureOrder"),
        "activeIndices": endpoint.get("activeIndices"),
        "endpointIntervalMargins": {
            "actualKrawczykQ": q,
            "scaledCompanionRadiusOverSimultaneousBudget": companion_ratio,
            "boundaryRowRadiusOverSimultaneousBudget": boundary_ratio,
            "endpointRadiusOverChartCapacity": chart_ratio,
            "allMarginsBelowOne": margins_ok,
        },
        "statuses": {
            "endpointCoefficientInputStatus": status(
                "synchronized endpoint coefficient input",
                endpoint_input_ok,
                (
                    "The synchronized 200-point endpoint certificate closes "
                    "center synchronization, analytic companion input, "
                    "analytic boundary input, and the Krawczyk input, with "
                    "all reported margin ratios below one."
                ),
            ),
            "persistentPlueckerChartStatus": status(
                "persistent Pluecker chart remains nonzero",
                pluecker_ok,
                (
                    "The exact endpoint map remains inside the certified "
                    "persistent Pluecker chart because the endpoint-radius "
                    "budget is below chart capacity."
                ),
            ),
            "endpointFullActiveRowRankStatus": status(
                "endpoint map has full active row rank",
                full_rank_ok,
                (
                    "Nonvanishing of the active Pluecker coordinate gives "
                    "full row rank for the exact active endpoint map."
                ),
            ),
            "endpointMapOntoActiveRowsStatus": status(
                "endpoint map is onto active rows",
                onto_ok,
                (
                    "A finite-dimensional linear map with full active row "
                    "rank is surjective onto the active row space."
                ),
            ),
            "boundedEndpointRightInverseStatus": status(
                "bounded endpoint right inverse",
                right_inverse_ok,
                (
                    "Surjectivity in finite dimension gives a bounded right "
                    "inverse.  The bound is controlled by the same certified "
                    "Pluecker separation margin."
                ),
            ),
        },
        "persistentPlueckerChartClosed": pluecker_ok,
        "endpointFullActiveRowRankClosed": full_rank_ok,
        "endpointMapOntoActiveRowsClosed": onto_ok,
        "boundedEndpointRightInverseClosed": right_inverse_ok,
        "proof": [
            "Import the synchronized 200-point coefficient/Krawczyk certificate.",
            "Check the companion, boundary, Krawczyk, and Pluecker margin ratios are below one.",
            "Conclude the exact endpoint map stays in the persistent Pluecker chart.",
            "Read off full active row rank, surjectivity onto active endpoint rows, and a bounded right inverse.",
        ],
        "remainingAnalyticGap": None if right_inverse_ok else "Endpoint full active row rank input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint full active row rank consequence theorem")
    print(f"  coefficient input: {endpoint_input_ok}")
    print(f"  Pluecker chart: {pluecker_ok}")
    print(f"  full active row rank: {full_rank_ok}")
    print(f"  bounded right inverse: {right_inverse_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
