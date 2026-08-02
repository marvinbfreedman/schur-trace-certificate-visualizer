#!/usr/bin/env python3
r"""Apply the continuous DAE inclusion theorem in the scaled crossing chart."""

from __future__ import annotations

import argparse
import json
import math
import sys
from decimal import Decimal, localcontext
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402
from flint import arb, ctx  # noqa: E402

from endpoint_flow_chebyshev_center import scaling_diag  # noqa: E402
from endpoint_restricted_dae_flow_center import base_shift  # noqa: E402
from endpoint_undivided_dae_krawczyk_enclosure import (  # noqa: E402
    arb_upper,
    ball,
    cheb_coefficients_ball,
    cheb_eval_and_derivative,
    point,
)
from endpoint_uniform_eigenrow_continuation_theorem import (  # noqa: E402
    arb_lower_decimal,
    arb_upper_decimal,
)
from global_trace_observability_gap import f  # noqa: E402


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": bool(closed),
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def interval_from_text(lower: str, upper: str) -> arb:
    lo = Decimal(lower)
    hi = Decimal(upper)
    return ball((lo + hi) / 2, (hi - lo) / 2)


def derivative_component(disk: dict, q: int, component: int) -> arb:
    row = disk["derivativeEnclosures"][q]["components"][component]
    return interval_from_text(
        row.get("lowerText", str(row["lower"])),
        row.get("upperText", str(row["upper"])),
    )


def adjoint_coefficients(disk: dict) -> list[arb]:
    out = []
    for m in range(9):
        value = arb(0)
        for k in range(m, 9):
            value += (
                point((-1) ** k * math.comb(k, m))
                * derivative_component(disk, k - m, k)
                / math.factorial(k)
            )
        out.append(value)
    return out


def matrix_inf(matrix: mp.matrix) -> mp.mpf:
    return max(
        mp.fsum(abs(matrix[i, j]) for j in range(matrix.cols))
        for i in range(matrix.rows)
    )


def interval_abs_upper(value: arb) -> mp.mpf:
    return mp.mpf(str(arb_upper_decimal(abs(value))))


def interval_abs_lower(value: arb) -> mp.mpf:
    if value.contains(0):
        return mp.mpf("0")
    lo = mp.mpf(str(arb_lower_decimal(value)))
    hi = mp.mpf(str(arb_upper_decimal(value)))
    return min(abs(lo), abs(hi))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cover-json", default="endpoint_eigenrow_between_node_crossing_cover.json"
    )
    parser.add_argument(
        "--graph-json", default="endpoint_undivided_dae_krawczyk_enclosure.json"
    )
    parser.add_argument(
        "--transmission-json", default="endpoint_dae_transmission_theorem.json"
    )
    parser.add_argument(
        "--inclusion-json", default="continuous_undivided_dae_residual_inclusion_theorem.json"
    )
    parser.add_argument(
        "--json-out", default="endpoint_crossing_continuous_residual_inclusion_certificate.json"
    )
    parser.add_argument("--bits", type=int, default=256)
    parser.add_argument("--scale-radius", default="0.2")
    parser.add_argument("--cell-subdivisions", type=int, default=32)
    args = parser.parse_args()

    ctx.prec = args.bits
    mp.mp.dps = max(80, int(args.bits * 0.32))
    cover = json.loads(Path(args.cover_json).read_text(encoding="utf-8"))
    graph = json.loads(Path(args.graph_json).read_text(encoding="utf-8"))
    transmission = json.loads(Path(args.transmission_json).read_text(encoding="utf-8"))
    inclusion = json.loads(Path(args.inclusion_json).read_text(encoding="utf-8"))
    crossing = next(row for row in graph["pieces"] if row["label"] == "crossing")
    left, right = map(mp.mpf, map(str, crossing["interval"]))
    nodes = [mp.mpf(str(value)) for value in crossing["lobattoNodes"]]
    scale_radius = mp.mpf(args.scale_radius)
    scale = mp.diag(
        [scale_radius**j / mp.factorial(j) for j in range(8)]
    )
    d7 = mp.diag([scale[j, j] for j in range(7)])
    bmat, bvec = base_shift()
    b_scaled = d7 * bmat * (d7 ** -1)
    bvec_scaled = d7 * bvec / scale[7, 7]
    b_norm = matrix_inf(b_scaled)
    bvec_norm = max(abs(bvec_scaled[i]) for i in range(7))

    raw_radius = mp.mpf(str(crossing["krawczykRadius"]))
    scaled_values = []
    scaled_radius = mp.mpf("0")
    for row in crossing["nodalCenters"]:
        scaled_row = []
        for j in range(7):
            ratio = scale[7, 7] / scale[j, j]
            scaled_row.append(ratio * mp.mpf(str(row[j])))
            scaled_radius = max(scaled_radius, abs(ratio) * raw_radius)
        scaled_values.append(scaled_row)
    coefficients = cheb_coefficients_ball(scaled_values, scaled_radius)
    midpoint = (left + right) / 2
    half_width = (right - left) / 2

    disks = cover["disks"]
    coefficient_rows = []
    alpha_lowers = []
    alpha_uppers = []
    c7_intervals = []
    for disk in disks:
        exact = adjoint_coefficients(disk)
        aprime = derivative_component(disk, 1, 8) / math.factorial(8)
        alpha_lowers.append(arb_lower_decimal(aprime))
        alpha_uppers.append(arb_upper_decimal(aprime))
        c7_intervals.append(exact[7])
        coefficient_rows.append((disk, exact))
    alpha = interval_from_text(
        format(min(alpha_lowers), ".80E"), format(max(alpha_uppers), ".80E")
    )
    alpha0 = interval_abs_lower(alpha)
    p0 = mp.inf
    for c7 in c7_intervals:
        ratio = c7 / alpha
        p0 = min(p0, mp.mpf(str(arb_lower_decimal(ratio))))

    residual_upper = mp.mpf("0")
    q_upper = mp.mpf("0")
    cell_rows = []
    for disk, exact in coefficient_rows:
        a = max(left, mp.mpf(str(disk["realCell"][0])))
        b = min(right, mp.mpf(str(disk["realCell"][1])))
        if a >= b:
            continue
        d_scaled = [
            point(scale[7, 7] / scale[k, k]) * exact[k] for k in range(7)
        ]
        cell_residual = mp.mpf("0")
        cell_q_upper = mp.mpf("0")
        for subdivision in range(args.cell_subdivisions):
            aa = a + (b - a) * subdivision / args.cell_subdivisions
            bb = a + (b - a) * (subdivision + 1) / args.cell_subdivisions
            sbox = ball((aa + bb) / 2, (bb - aa) / 2)
            xbox = (sbox - point(midpoint)) / point(half_width)
            q, qprime = cheb_eval_and_derivative(coefficients, xbox, half_width)
            local_q_upper = max(interval_abs_upper(value) for value in q)
            cell_q_upper = max(cell_q_upper, local_q_upper)
            qb_scalar = sum(
                (q[j] * point(bvec_scaled[j]) for j in range(7)), arb(0)
            )
            for k in range(7):
                qb = sum(
                    (q[j] * point(b_scaled[j, k]) for j in range(7)), arb(0)
                )
                residual = (
                    exact[8] * (qprime[k] + qb + qb_scalar * q[k])
                    + exact[7] * q[k]
                    + d_scaled[k]
                )
                cell_residual = max(cell_residual, interval_abs_upper(residual))
        q_upper = max(q_upper, cell_q_upper)
        residual_upper = max(residual_upper, cell_residual)
        cell_rows.append(
            {
                "interval": [f(a), f(b)],
                "qScaledUpper": f(cell_q_upper),
                "exactResidualUpper": f(cell_residual),
            }
        )

    root_left, root_right = map(
        mp.mpf, map(str, transmission["domain"]["certifiedZeroBracket"])
    )
    root_mid = (root_left + root_right) / 2
    h = max(root_mid - left, right - root_mid) + (root_right - root_left) / 2
    eta = residual_upper / (alpha0 * p0) if alpha0 > 0 and p0 > 0 else mp.inf
    kappa = h * (b_norm + 2 * bvec_norm * q_upper) / p0 if p0 > 0 else mp.inf
    nu = 2 * h * bvec_norm / p0 if p0 > 0 else mp.inf
    discriminant = (1 - kappa) ** 2 - 4 * nu * eta
    radius = mp.inf
    if discriminant > 0 and kappa < 1:
        small_root = 2 * eta / ((1 - kappa) + mp.sqrt(discriminant))
        radius = 2 * small_root
    invariant = (
        mp.isfinite(radius)
        and eta + kappa * radius + nu * radius**2 < radius
    )
    contraction = kappa + 2 * nu * radius if mp.isfinite(radius) else mp.inf
    closed = (
        cover.get("betweenNodeEigenrowDerivativeCoverClosed", False)
        and inclusion.get("continuousResidualInclusionTheoremClosed", False)
        and alpha0 > 0
        and p0 > 0
        and invariant
        and contraction < 1
    )
    data = {
        "schemaVersion": 1,
        "theoremName": "exact continuous residual inclusion in the scaled endpoint crossing chart",
        "proofClass": "interval/ball certificate plus analytic residual inclusion",
        "statement": (
            "The exact undivided endpoint graph DAE on the crossing interval has a unique "
            "bounded solution in the displayed scaled sup-norm tube about the certified "
            "Chebyshev graph whenever the final status is closed."
        ),
        "domain": {
            "interval": [f(left), f(right)],
            "stateScaling": "Taylor diagonal scale radius^j/j!",
            "scaleRadius": f(scale_radius),
            "cellSubdivisions": args.cell_subdivisions,
            "crossingRootBracket": [f(root_left), f(root_right)],
        },
        "dependencies": [
            args.cover_json,
            args.graph_json,
            args.transmission_json,
            args.inclusion_json,
        ],
        "scaledConstants": {
            "BInfNorm": f(b_norm),
            "bInfNorm": f(bvec_norm),
            "qCandidateSup": f(q_upper),
            "graphNodalRadius": f(scaled_radius),
            "crossingHalfLengthFromRoot": f(h),
            "alphaInterval": [
                format(arb_lower_decimal(alpha), ".80E"),
                format(arb_upper_decimal(alpha), ".80E"),
            ],
            "alphaAbsLower": f(alpha0),
            "pLower": f(p0),
            "exactUndividedResidualUpper": f(residual_upper),
        },
        "radiiPolynomial": {
            "eta": f(eta),
            "kappa": f(kappa),
            "nu": f(nu),
            "discriminant": f(discriminant),
            "selectedRadius": f(radius),
            "contractionBound": f(contraction),
            "invariantBall": bool(invariant),
        },
        "cells": cell_rows,
        "continuousCrossingResidualInclusionClosed": closed,
        "continuousCrossingResidualInclusionStatus": status(
            "exact continuous crossing residual-to-solution inclusion",
            closed,
            (
                "The exact coefficient cover and scaled singular radii polynomial close."
                if closed
                else "The displayed exact residual or scaled singular radii inequality does not close."
            ),
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"alpha0={mp.nstr(alpha0, 10)} p0={mp.nstr(p0, 10)}")
    print(f"Q={mp.nstr(q_upper, 10)} residual={mp.nstr(residual_upper, 10)}")
    print(
        f"eta={mp.nstr(eta, 10)} kappa={mp.nstr(kappa, 10)} "
        f"nu={mp.nstr(nu, 10)} contraction={mp.nstr(contraction, 10)}"
    )
    print(f"continuous crossing inclusion closed: {closed}")
    print(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
