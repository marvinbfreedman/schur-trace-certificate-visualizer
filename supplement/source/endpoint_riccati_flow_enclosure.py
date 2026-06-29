#!/usr/bin/env python3
r"""Riccati/Pluecker chart enclosure for the endpoint full-rank shortcut.

This is the ODE layer under the chart-ball certificate.  For the row flow

    M'(s) = -M(s) A(s),

the Pluecker vector p_ij = det M[:,(i,j)] satisfies a linear exterior-square
ODE

    p'(s) = B(s) p(s).

In the persistent chart p_01 != 0, the ratios

    z_J = p_J / p_01

satisfy a Riccati equation.  This script constructs that exact finite-dimensional
Riccati vector field from the same confluent trace derivatives used by the
endpoint-flow center, then combines it with the existing chart-ball certificate.

The output is intentionally precise about status:

* the exterior/Riccati equation and chart-ball implication are closed;
* the final exact-flow enclosure is closed only after a validated residual or
  interval integration proves that the exact normalized Pluecker vector lies in
  the certified projective ball.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import companion_matrix, status  # noqa: E402
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    aligned_derivative_cache,
    node_key,
    parse_ints,
)
from endpoint_flow_chebyshev_center import cheb_lobatto_nodes  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def pairs(dim: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(dim) for j in range(i + 1, dim)]


def pair_index_map(dim: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(pairs(dim))}


def antisym_pair(a: int, b: int, index: dict[tuple[int, int], int]):
    if a == b:
        return None, 0
    if a < b:
        return index[(a, b)], 1
    return index[(b, a)], -1


def exterior_generator_for_row_flow(amat: mp.matrix) -> mp.matrix:
    """Return B for p'=B p when row vectors satisfy r'=-r A."""
    dim = amat.rows
    ps = pairs(dim)
    index = pair_index_map(dim)
    out = mp.matrix(len(ps), len(ps))
    for row_idx, (i, j) in enumerate(ps):
        for k in range(dim):
            # Contribution -A[k,i] p_{k,j}.
            col_idx, sign = antisym_pair(k, j, index)
            if col_idx is not None:
                out[row_idx, col_idx] += -amat[k, i] * sign
            # Contribution -A[k,j] p_{i,k}.
            col_idx, sign = antisym_pair(i, k, index)
            if col_idx is not None:
                out[row_idx, col_idx] += -amat[k, j] * sign
    return out


def plucker_vector(mat: list[list[float]]) -> list[mp.mpf]:
    out = []
    dim = len(mat[0])
    for i, j in pairs(dim):
        out.append(
            mp.mpf(str(mat[0][i])) * mp.mpf(str(mat[1][j]))
            - mp.mpf(str(mat[0][j])) * mp.mpf(str(mat[1][i]))
        )
    return out


def normalize(vec: list[mp.mpf]) -> list[mp.mpf]:
    norm = mp.sqrt(mp.fsum(value * value for value in vec))
    return [value / norm for value in vec]


def orient_chart(vec: list[mp.mpf], chart_idx: int) -> list[mp.mpf]:
    if vec[chart_idx] < 0:
        return [-value for value in vec]
    return vec


def ratios_from_normalized(vec: list[mp.mpf], chart_idx: int) -> list[mp.mpf]:
    pivot = vec[chart_idx]
    return [value / pivot for value in vec]


def riccati_rhs_from_ratios(bmat: mp.matrix, ratios: list[mp.mpf], chart_idx: int) -> list[mp.mpf]:
    scalar = mp.fsum(bmat[chart_idx, k] * ratios[k] for k in range(bmat.cols))
    out = []
    for j in range(bmat.rows):
        value = mp.fsum(bmat[j, k] * ratios[k] for k in range(bmat.cols))
        out.append(value - ratios[j] * scalar)
    return out


def matrix_abs_max(mats: list[mp.matrix]) -> mp.matrix:
    rows = mats[0].rows
    cols = mats[0].cols
    out = mp.matrix(rows, cols)
    for mat in mats:
        for i in range(rows):
            for j in range(cols):
                out[i, j] = max(out[i, j], abs(mat[i, j]))
    return out


def row_sum_norm(mat: mp.matrix) -> mp.mpf:
    return max(mp.fsum(abs(mat[i, j]) for j in range(mat.cols)) for i in range(mat.rows))


def lipschitz_bound(babs: mp.matrix, zabs: list[mp.mpf], chart_idx: int) -> mp.mpf:
    """Infinity-norm Lipschitz bound for the chart Riccati vector field."""
    scalar_abs = mp.fsum(babs[chart_idx, k] * zabs[k] for k in range(babs.cols))
    max_row = mp.mpf("0")
    for j in range(babs.rows):
        total = mp.mpf("0")
        for ell in range(babs.cols):
            value = babs[j, ell] + zabs[j] * babs[chart_idx, ell]
            if j == ell:
                value += scalar_abs
            total += value
        max_row = max(max_row, total)
    return max_row


def rhs_bound(babs: mp.matrix, zabs: list[mp.mpf], chart_idx: int) -> mp.mpf:
    scalar_abs = mp.fsum(babs[chart_idx, k] * zabs[k] for k in range(babs.cols))
    max_rhs = mp.mpf("0")
    for j in range(babs.rows):
        linear = mp.fsum(babs[j, k] * zabs[k] for k in range(babs.cols))
        max_rhs = max(max_rhs, linear + zabs[j] * scalar_abs)
    return max_rhs


def chart_ratio_bounds(chart_ball: dict) -> list[mp.mpf]:
    bounds = [
        mp.mpf(str(row["absRatioUpper"]))
        for row in chart_ball["certifiedProjectiveBall"]["chartRatioBounds"]
    ]
    return bounds


def build_nodes(args) -> list[mp.mpf]:
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    nodes = []
    for order in parse_ints(args.orders):
        nodes.extend(cheb_lobatto_nodes(order, left, s0))
        nodes.extend(cheb_lobatto_nodes(order, s0, right))
    unique = sorted({node_key(node): node for node in nodes}.values())
    return unique


def derivative_args(args):
    return SimpleNamespace(**vars(args))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grassmann-json", default="endpoint_grassmann_flow_center.json")
    parser.add_argument("--chart-ball-json", default="endpoint_grassmann_chart_ball_certificate.json")
    parser.add_argument("--row-flow-json", default="endpoint_adjoint_row_flow_center_7_9_11.json")
    parser.add_argument("--orders", default="7 9 11")
    parser.add_argument("--coefficient-inflation", default="1e-10")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=11)
    parser.add_argument("--active-tol", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="endpoint_riccati_flow_enclosure.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    grassmann = load_json(args.grassmann_json)
    chart_ball = load_json(args.chart_ball_json)
    row_flow = load_json(args.row_flow_json)
    chart_idx = int(grassmann["persistentChart"]["coordinateIndex"])
    chart_cols = grassmann["persistentChart"]["columns"]

    latest_center = grassmann["grassmannCenters"][-1]
    center_vec = [
        mp.mpf(str(value))
        for value in latest_center["orientedNormalizedVector"]
    ]
    center_vec = orient_chart(center_vec, chart_idx)
    center_ratios = ratios_from_normalized(center_vec, chart_idx)
    ratio_bounds = chart_ratio_bounds(chart_ball)

    nodes = build_nodes(args)
    max_q = max(args.max_trace_q, args.jet_order - 1)
    print(
        f"Building exact exterior/Riccati generators on {len(nodes)} nodes "
        f"max_q={max_q}",
        flush=True,
    )
    cache, derivative_rows, min_gap, min_cos = aligned_derivative_cache(
        derivative_args(args),
        nodes,
        max_q,
    )
    b_mats = []
    rhs_rows = []
    for node in nodes:
        derivs = cache[node_key(node)]
        amat = companion_matrix(derivs, args.jet_order - 1)
        bmat = exterior_generator_for_row_flow(amat)
        b_mats.append(bmat)
        rhs = riccati_rhs_from_ratios(bmat, center_ratios, chart_idx)
        rhs_rows.append(
            {
                "s": f(node),
                "riccatiRhsInfAtCenter": f(max(abs(value) for value in rhs)),
                "exteriorGeneratorInfNorm": f(row_sum_norm(bmat)),
            }
        )

    babs = matrix_abs_max(b_mats)
    inflation = mp.mpf(args.coefficient_inflation)
    babs_inflated = (1 + inflation) * babs
    lipschitz = lipschitz_bound(babs_inflated, ratio_bounds, chart_idx)
    rhs = rhs_bound(babs_inflated, ratio_bounds, chart_idx)
    interval_length = mp.mpf(args.constraint_max) - mp.mpf(args.constraint_min)
    cert_radius = mp.mpf(str(chart_ball["certifiedProjectiveRadius"]))
    threshold_capacity = mp.mpf(str(chart_ball["radiusCapacityForThreshold"]))
    latest_projective = mp.mpf(str(chart_ball["latestProjectiveDistance"]))
    residual_budget = cert_radius / interval_length
    if lipschitz > 0:
        gronwall_residual_budget = cert_radius * lipschitz / (mp.e ** (lipschitz * interval_length) - 1)
    else:
        gronwall_residual_budget = residual_budget

    row_flow_closed = bool(
        latest_projective * 4 < cert_radius
        and cert_radius < threshold_capacity
    )
    data = {
        "theoremName": "endpoint Riccati flow enclosure",
        "sourceGrassmannJson": args.grassmann_json,
        "sourceChartBallJson": args.chart_ball_json,
        "sourceRowFlowJson": args.row_flow_json,
        "chartCoordinate": chart_cols,
        "chartIndex": chart_idx,
        "orders": parse_ints(args.orders),
        "nodeCount": len(nodes),
        "coefficientInflation": f(inflation),
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "centerOrder": latest_center["chebyshevOrder"],
        "centerPivotAbs": latest_center["chartMinor"]["absNormalizedCoordinate"],
        "certifiedProjectiveRadius": f(cert_radius),
        "radiusCapacityForThreshold": f(threshold_capacity),
        "latestProjectiveDistance": f(latest_projective),
        "fourTimesLatestProjectiveDistance": f(4 * latest_projective),
        "ratioBoundsFromChartBall": [f(value) for value in ratio_bounds],
        "exteriorGeneratorInfNormBound": f(row_sum_norm(babs_inflated)),
        "riccatiRhsInfBoundOnChartTube": f(rhs),
        "riccatiLipschitzInfBoundOnChartTube": f(lipschitz),
        "intervalLength": f(interval_length),
        "flatResidualBudgetForRadius": f(residual_budget),
        "gronwallResidualBudgetForRadius": f(gronwall_residual_budget),
        "gronwallResidualBudgetForRadiusText": mp.nstr(gronwall_residual_budget, 16),
        "generatorSamples": rhs_rows,
        "exteriorEquationStatus": status(
            "exact exterior-square equation",
            True,
            (
                "For row flow M'=-MA, the Pluecker vector satisfies p'=B(s)p "
                "with B constructed by the exterior-square formula in this "
                "certificate."
            ),
        ),
        "riccatiChartEquationStatus": status(
            "exact Riccati chart equation",
            True,
            (
                "In the chart p_01 != 0, ratios z_J=p_J/p_01 satisfy the "
                "displayed Riccati equation z'=Bp/p_01-z(Bp)_01/p_01."
            ),
        ),
        "chartTubeNoncollapseStatus": status(
            "chart tube noncollapse",
            bool(cert_radius < threshold_capacity),
            (
                "The certified projective tube lies within the chart-ball "
                "radius needed for |p_hat_01| >= 1/2."
            ),
        ),
        "empiricalRefinementInsideTubeStatus": status(
            "empirical refinement inside chart tube",
            row_flow_closed,
            (
                "The latest observed 9->11 projective movement, inflated by "
                "a factor four, is inside the certified projective radius. "
                "This is evidence, not the exact residual proof."
            ),
        ),
        "validatedExactFlowEnclosureStatus": status(
            "validated exact endpoint Riccati flow enclosure",
            False,
            (
                "Open.  The ODE and chart tube are now explicit.  To close "
                "this, prove a residual or interval-integration bound below "
                "the listed projective radius budget."
            ),
        ),
        "conclusion": (
            "The exact exterior/Riccati equation has been constructed in the "
            "persistent p_01 chart.  The chart tube has large noncollapse "
            "margin; the remaining proof obligation is a validated residual "
            "or interval integration bound for the exact endpoint BVP flow."
        ),
        "nextTarget": (
            "Use a Chebyshev/Krawczyk or interval Taylor solve for the Riccati "
            "BVP and prove the projective enclosure radius is below the "
            "certified chart-ball radius."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Endpoint Riccati flow enclosure")
    print(f"  chart: {chart_cols}")
    print(f"  nodes: {len(nodes)}")
    print(f"  center |p_hat_chart|: {latest_center['chartMinor']['absNormalizedCoordinate']:.12e}")
    print(f"  certified projective radius: {fmt(cert_radius, 12)}")
    print(f"  4x latest projective distance: {fmt(4 * latest_projective, 12)}")
    print(f"  generator inf-norm bound: {fmt(row_sum_norm(babs_inflated), 12)}")
    print(f"  Riccati RHS bound on tube: {fmt(rhs, 12)}")
    print(f"  Riccati Lipschitz bound on tube: {fmt(lipschitz, 12)}")
    print(f"  flat residual budget: {fmt(residual_budget, 12)}")
    print(f"  Gronwall residual budget: {fmt(gronwall_residual_budget, 12)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
