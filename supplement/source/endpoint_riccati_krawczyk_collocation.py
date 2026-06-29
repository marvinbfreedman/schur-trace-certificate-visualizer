#!/usr/bin/env python3
r"""Krawczyk/interval-collocation target for the endpoint Riccati chart.

The endpoint Pluecker chart point is obtained from two row-flow collocation
problems:

    r'(s) = -r(s)A(s),

one transporting the left endpoint rows to s0 and one transporting the right
endpoint rows to s0.  The endpoint map is their difference, and the persistent
Riccati chart is p_01 != 0.

This script builds the finite Chebyshev collocation systems, computes the
Krawczyk constants for coefficient/boundary perturbations, and converts the
resulting endpoint-entry radius into the already certified Pluecker chart
ball.  It is the interval-collocation proof skeleton: once the coefficient and
boundary rows are enclosed below the reported radii, the chart noncollapse
follows.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    companion_matrix,
    status,
)
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    active_endpoint_covectors,
    aligned_derivative_cache,
    nearest_derivs,
    node_key,
    scaling_diag,
)
from endpoint_flow_chebyshev_center import (  # noqa: E402
    cheb_lobatto_diff_matrix,
    cheb_lobatto_nodes,
)
from endpoint_grassmann_chart_ball_certificate import raw_entry_ratio_lower  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def matrix_inf_norm(mat: mp.matrix) -> mp.mpf:
    return max(mp.fsum(abs(mat[i, j]) for j in range(mat.cols)) for i in range(mat.rows))


def matrix_max_abs(mat: mp.matrix) -> mp.mpf:
    return max(abs(mat[i, j]) for i in range(mat.rows) for j in range(mat.cols))


def build_linear_system(
    nodes: list[mp.mpf],
    derivs,
    boundary_index: int,
    boundary_rows: mp.matrix,
    scaling: str,
):
    order = len(nodes)
    dim = boundary_rows.cols
    left = min(nodes)
    right = max(nodes)
    dmat = cheb_lobatto_diff_matrix(order, left, right)
    scale = scaling_diag(nodes, scaling, dim)
    scale_inv = scale ** -1
    amats = [
        scale * companion_matrix(derivs[j], dim) * scale_inv
        for j in range(order)
    ]
    boundary_scaled = boundary_rows * scale_inv
    size = order * dim
    system = mp.matrix(size, size)
    rhs = mp.matrix(size, boundary_rows.rows)

    def idx(node: int, comp: int) -> int:
        return node * dim + comp

    for j in range(order):
        for m in range(dim):
            row = idx(j, m)
            if j == boundary_index:
                system[row, idx(j, m)] = 1
                for r in range(boundary_rows.rows):
                    rhs[row, r] = boundary_scaled[r, m]
                continue
            for k in range(order):
                system[row, idx(k, m)] += dmat[j, k]
            amat = amats[j]
            for ell in range(dim):
                system[row, idx(j, ell)] += amat[ell, m]
    return system, rhs, scale


def solve_system(system: mp.matrix, rhs: mp.matrix) -> mp.matrix:
    out = mp.matrix(system.rows, rhs.cols)
    for col in range(rhs.cols):
        sol = mp.lu_solve(system, rhs[:, col])
        for row in range(system.rows):
            out[row, col] = sol[row]
    return out


def target_rows(solution: mp.matrix, target_index: int, dim: int, scale: mp.matrix) -> mp.matrix:
    out = mp.matrix(solution.cols, dim)
    for r in range(solution.cols):
        for m in range(dim):
            out[r, m] = solution[target_index * dim + m, r] * scale[m, m]
    return out


def krawczyk_bound(inv_norm, system_delta, rhs_delta, solution_norm):
    q = inv_norm * system_delta
    if q >= 1:
        return mp.inf, q
    err = inv_norm * (rhs_delta + system_delta * solution_norm) / (1 - q)
    return err, q


def endpoint_entry_radius(left: dict, right: dict, coeff_radius: mp.mpf, boundary_radius: mp.mpf):
    dim = left["dim"]
    delta_s = dim * coeff_radius
    left_rhs_delta = boundary_radius * left["scaleInvInf"]
    right_rhs_delta = boundary_radius * right["scaleInvInf"]
    left_err, left_q = krawczyk_bound(
        left["inverseInfNorm"],
        delta_s,
        left_rhs_delta,
        left["solutionInfNorm"],
    )
    right_err, right_q = krawczyk_bound(
        right["inverseInfNorm"],
        delta_s,
        right_rhs_delta,
        right["solutionInfNorm"],
    )
    scale_out = max(left["scaleInf"], right["scaleInf"])
    return scale_out * (left_err + right_err), max(left_q, right_q)


def coeff_radius_capacity(left, right, boundary_radius, target_entry_radius):
    hi = mp.mpf("1")
    while True:
        radius, q = endpoint_entry_radius(left, right, hi, boundary_radius)
        if q >= 1 or radius > target_entry_radius:
            break
        hi *= 2
        if hi > mp.mpf("1e50"):
            break
    lo = mp.mpf("0")
    for _ in range(100):
        mid = (lo + hi) / 2
        radius, q = endpoint_entry_radius(left, right, mid, boundary_radius)
        if q < 1 and radius <= target_entry_radius:
            lo = mid
        else:
            hi = mid
    radius, q = endpoint_entry_radius(left, right, lo, boundary_radius)
    return lo, radius, q


def boundary_radius_capacity(left, right, coeff_radius, target_entry_radius):
    hi = mp.mpf("1")
    while True:
        radius, q = endpoint_entry_radius(left, right, coeff_radius, hi)
        if q >= 1 or radius > target_entry_radius:
            break
        hi *= 2
        if hi > mp.mpf("1e50"):
            break
    lo = mp.mpf("0")
    for _ in range(100):
        mid = (lo + hi) / 2
        radius, q = endpoint_entry_radius(left, right, coeff_radius, mid)
        if q < 1 and radius <= target_entry_radius:
            lo = mid
        else:
            hi = mid
    radius, q = endpoint_entry_radius(left, right, coeff_radius, lo)
    return lo, radius, q


def simultaneous_fraction_capacity(left, right, coeff_cap, boundary_cap, target_entry_radius):
    hi = mp.mpf("1")
    lo = mp.mpf("0")
    for _ in range(100):
        mid = (lo + hi) / 2
        radius, q = endpoint_entry_radius(
            left,
            right,
            mid * coeff_cap,
            mid * boundary_cap,
        )
        if q < 1 and radius <= target_entry_radius:
            lo = mid
        else:
            hi = mid
    radius, q = endpoint_entry_radius(
        left,
        right,
        lo * coeff_cap,
        lo * boundary_cap,
    )
    return lo, radius, q


def side_certificate(system, rhs, solution, scale, target_index, dim, label):
    inv = system ** -1
    residual = system * solution - rhs
    return {
        "label": label,
        "dim": dim,
        "systemSize": system.rows,
        "systemInfNorm": matrix_inf_norm(system),
        "inverseInfNorm": matrix_inf_norm(inv),
        "conditionInfEstimate": matrix_inf_norm(system) * matrix_inf_norm(inv),
        "solutionInfNorm": matrix_max_abs(solution),
        "residualInfNorm": matrix_max_abs(residual),
        "scaleInf": max(abs(scale[i, i]) for i in range(scale.rows)),
        "scaleInvInf": max(abs(1 / scale[i, i]) for i in range(scale.rows)),
        "targetIndex": target_index,
    }


def floatside(row: dict) -> dict:
    return {
        key: f(value) if isinstance(value, mp.mpf) else value
        for key, value in row.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--chart-ball-json", default="endpoint_grassmann_chart_ball_certificate.json")
    parser.add_argument("--row-flow-json", default="endpoint_adjoint_row_flow_center_7_9_11.json")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
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
    parser.add_argument("--json-out", default="endpoint_riccati_krawczyk_collocation.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    chart_ball = load_json(args.chart_ball_json)
    row_flow = load_json(args.row_flow_json)
    raw_center = next(
        row for row in row_flow["centers"]
        if row["chebyshevOrder"] == args.order
    )
    chart_cols = chart_ball["chartCoordinate"]
    raw_capacity = mp.mpf(str(chart_ball["rawEntryBallCapacity"]["uniformAbsoluteEntryRadius"]))

    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    dim = args.jet_order - 1
    max_q = max(args.max_trace_q, dim)
    left_nodes = cheb_lobatto_nodes(args.order, left, s0)
    right_nodes = cheb_lobatto_nodes(args.order, s0, right)
    all_nodes = left_nodes + right_nodes
    print(f"Building Krawczyk collocation systems order={args.order}", flush=True)
    cache, derivative_rows, min_gap, min_cos = aligned_derivative_cache(
        SimpleNamespace(**vars(args)),
        all_nodes,
        max_q,
    )
    left_derivs = [cache[node_key(node)] for node in left_nodes]
    right_derivs = [cache[node_key(node)] for node in right_nodes]
    e_s0 = nearest_derivs(cache, s0)
    e_left = nearest_derivs(cache, left)
    e_right = nearest_derivs(cache, right)
    base, active_modes, active_idx = active_setup(args, e_s0)
    c_left = active_endpoint_covectors(base["polys"], active_modes, left, e_left, args.jet_order)
    c_right = active_endpoint_covectors(base["polys"], active_modes, right, e_right, args.jet_order)

    left_system, left_rhs, left_scale = build_linear_system(
        left_nodes,
        left_derivs,
        args.order - 1,
        c_left,
        args.state_scaling,
    )
    right_system, right_rhs, right_scale = build_linear_system(
        right_nodes,
        right_derivs,
        0,
        c_right,
        args.state_scaling,
    )
    left_solution = solve_system(left_system, left_rhs)
    right_solution = solve_system(right_system, right_rhs)
    left_target = target_rows(left_solution, 0, dim, left_scale)
    right_target = target_rows(right_solution, args.order - 1, dim, right_scale)
    endpoint = right_target - left_target
    endpoint_delta = max(
        abs(endpoint[i, j] - mp.mpf(str(raw_center["endpointMap"][i][j])))
        for i in range(endpoint.rows)
        for j in range(endpoint.cols)
    )
    left_cert = side_certificate(
        left_system,
        left_rhs,
        left_solution,
        left_scale,
        0,
        dim,
        "left-to-s0 row flow",
    )
    right_cert = side_certificate(
        right_system,
        right_rhs,
        right_solution,
        right_scale,
        args.order - 1,
        dim,
        "right-to-s0 row flow",
    )
    coeff_cap, coeff_entry_radius, coeff_q = coeff_radius_capacity(
        left_cert,
        right_cert,
        mp.mpf("0"),
        raw_capacity,
    )
    boundary_cap, boundary_entry_radius, boundary_q = boundary_radius_capacity(
        left_cert,
        right_cert,
        mp.mpf("0"),
        raw_capacity,
    )
    half_coeff_radius, half_coeff_q = endpoint_entry_radius(
        left_cert,
        right_cert,
        coeff_cap / 2,
        boundary_cap / 2,
    )
    half_ratio_lower = raw_entry_ratio_lower(
        raw_center["endpointMap"],
        chart_cols,
        half_coeff_radius,
    )
    simultaneous_factor, simultaneous_radius, simultaneous_q = simultaneous_fraction_capacity(
        left_cert,
        right_cert,
        coeff_cap,
        boundary_cap,
        raw_capacity,
    )
    simultaneous_ratio_lower = raw_entry_ratio_lower(
        raw_center["endpointMap"],
        chart_cols,
        simultaneous_radius,
    )
    data = {
        "theoremName": "endpoint Riccati Krawczyk collocation certificate",
        "order": args.order,
        "stateScaling": args.state_scaling,
        "chartCoordinate": chart_cols,
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "rawEndpointEntryCapacityForChartThreshold": f(raw_capacity),
        "recomputedEndpointMaxDifference": f(endpoint_delta),
        "leftSystem": floatside(left_cert),
        "rightSystem": floatside(right_cert),
        "coefficientRadiusCapacityScaledCompanion": f(coeff_cap),
        "endpointEntryRadiusAtCoefficientCapacity": f(coeff_entry_radius),
        "krawczykQAtCoefficientCapacity": f(coeff_q),
        "boundaryRadiusCapacityUnscaledRows": f(boundary_cap),
        "endpointEntryRadiusAtBoundaryCapacity": f(boundary_entry_radius),
        "krawczykQAtBoundaryCapacity": f(boundary_q),
        "mixedHalfCapacityEndpointEntryRadius": f(half_coeff_radius),
        "mixedHalfCapacityKrawczykQ": f(half_coeff_q),
        "mixedHalfCapacityPlueckerRatioLower": f(half_ratio_lower),
        "simultaneousCapacityFraction": f(simultaneous_factor),
        "simultaneousCoefficientRadius": f(simultaneous_factor * coeff_cap),
        "simultaneousBoundaryRadius": f(simultaneous_factor * boundary_cap),
        "simultaneousEndpointEntryRadius": f(simultaneous_radius),
        "simultaneousKrawczykQ": f(simultaneous_q),
        "simultaneousPlueckerRatioLower": f(simultaneous_ratio_lower),
        "linearCollocationKrawczykStatus": status(
            "linear row-flow Krawczyk theorem",
            bool(coeff_q < 1 and boundary_q < 1),
            (
                "The finite collocation systems have explicit inverse norms. "
                "If scaled companion-matrix and endpoint-row interval radii "
                "are below the reported capacities, the endpoint map remains "
                "inside the persistent Pluecker chart threshold."
            ),
        ),
        "exactCoefficientEnclosureStatus": status(
            "exact coefficient and boundary interval enclosure",
            False,
            (
                "Open.  The Krawczyk radii are now quantified, but the exact "
                "confluent trace-derivative coefficients and boundary rows "
                "still need interval/ball enclosures below these radii."
            ),
        ),
        "conclusion": (
            "The Krawczyk/collocation proof has been reduced to finite "
            "coefficient and boundary-row interval radii.  This avoids the "
            "useless raw Riccati Gronwall bound, but exact coefficient balls "
            "are still required to close the theorem."
        ),
        "nextTarget": (
            "Build interval balls for the scaled companion matrices and active "
            "endpoint boundary rows, then compare their radii with the "
            "Krawczyk capacities reported here."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint Riccati Krawczyk collocation certificate")
    print(f"  order: {args.order}")
    print(f"  recomputed endpoint diff: {fmt(endpoint_delta, 12)}")
    print(f"  raw endpoint entry capacity: {fmt(raw_capacity, 12)}")
    print(f"  left cond_inf: {fmt(left_cert['conditionInfEstimate'], 12)}")
    print(f"  right cond_inf: {fmt(right_cert['conditionInfEstimate'], 12)}")
    print(f"  coefficient radius capacity: {fmt(coeff_cap, 12)}")
    print(f"  boundary row radius capacity: {fmt(boundary_cap, 12)}")
    print(f"  mixed half-capacity endpoint radius: {fmt(half_coeff_radius, 12)}")
    print(f"  mixed half-capacity |p_hat_chart| lower: {fmt(half_ratio_lower, 12)}")
    print(f"  simultaneous capacity fraction: {fmt(simultaneous_factor, 12)}")
    print(f"  simultaneous |p_hat_chart| lower: {fmt(simultaneous_ratio_lower, 12)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
