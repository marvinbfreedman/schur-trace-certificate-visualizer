#!/usr/bin/env python3
r"""Synchronized endpoint coefficient/Krawczyk input certificate.

This is the bookkeeping layer after the segment Bernstein and eigenrow
Taylor certificates.  It rebuilds the finite row-flow Krawczyk inputs at the
same Gauss-Legendre center controlled by the deterministic segment theorem
and compares analytic input radii against the Krawczyk capacities:

* scaled companion matrices, using the rational companion-map derivative;
* active endpoint boundary rows, using the linear Green-concomitant formula.

The active source frame is fixed at the synchronized center.  This is the
right finite-coordinate statement for the endpoint full-rank shortcut.
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
    adjoint_coefficients,
    companion_matrix,
    status,
)
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    active_endpoint_covectors,
    node_key,
    scaling_diag,
)
from endpoint_confluent_trace_tail_certificate import collocation_nodes  # noqa: E402
from endpoint_eigenrow_interval_propagation import build_center as eigenrow_center  # noqa: E402
from endpoint_flow_chebyshev_center import cheb_lobatto_nodes, exact_derivative_at  # noqa: E402
from endpoint_grassmann_chart_ball_certificate import (  # noqa: E402
    raw_entry_radius_capacity,
    raw_entry_ratio_lower,
)
from endpoint_riccati_krawczyk_collocation import (  # noqa: E402
    boundary_radius_capacity,
    build_linear_system,
    coeff_radius_capacity,
    endpoint_entry_radius,
    floatside,
    side_certificate,
    simultaneous_fraction_capacity,
    solve_system,
    target_rows,
)
from global_trace_observability_gap import f, fmt  # noqa: E402
from quotient_factorization_mp import poly_derivative_value  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def dot(a, b) -> mp.mpf:
    return mp.fsum(a[i] * b[i] for i in range(len(a)))


def flip_derivs(e_derivs):
    return [[-value for value in row] for row in e_derivs]


def matrix_to_lists(mat: mp.matrix) -> list[list[float]]:
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def matrix_max_abs(mat: mp.matrix) -> mp.mpf:
    return max(abs(mat[i, j]) for i in range(mat.rows) for j in range(mat.cols))


def max_abs_matrix(mats: list[mp.matrix]) -> mp.mpf:
    return max(matrix_max_abs(mat) for mat in mats) if mats else mp.mpf("0")


def derivative_cache(args, nodes: list[mp.mpf], max_q: int):
    unique = sorted({node_key(node): node for node in nodes}.values())
    cache = {}
    previous = None
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    rows = []
    hits = 0
    for idx, node in enumerate(unique):
        _mats, vals, e_derivs, _lam_derivs, cache_hit = eigenrow_center(args, node, max_q)
        hits += 1 if cache_hit else 0
        if previous is not None:
            cos = dot(previous, e_derivs[0])
            if cos < 0:
                e_derivs = flip_derivs(e_derivs)
                cos = -cos
            min_cos = min(min_cos, abs(cos))
        previous = list(e_derivs[0])
        gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
        min_gap = min(min_gap, gap)
        cache[node_key(node)] = e_derivs
        rows.append(
            {
                "index": idx,
                "s": f(node),
                "gap": f(gap),
                "cacheHit": cache_hit,
            }
        )
    return cache, rows, min_gap, min_cos, hits


def nearest_derivs(cache: dict[str, list[list[mp.mpf]]], value: mp.mpf):
    key = min(cache, key=lambda candidate: abs(mp.mpf(candidate) - value))
    return cache[key]


def companion_entry_radius(e_derivs, delta: mp.mpf, dim: int, scale: mp.matrix) -> mp.mpf:
    coeffs = adjoint_coefficients(e_derivs, dim)
    lead = abs(coeffs[dim])
    coeff_delta = []
    for m in range(dim + 1):
        coeff_delta.append(
            delta
            * mp.fsum(mp.binomial(k, m) / mp.factorial(k) for k in range(m, dim + 1))
        )
    lead_delta = coeff_delta[dim]
    lead_low = lead - lead_delta
    if lead_low <= 0:
        return mp.inf
    out = mp.mpf("0")
    for j in range(dim):
        ratio_delta = (
            coeff_delta[j] * lead + abs(coeffs[j]) * lead_delta
        ) / (lead_low * lead)
        scaled = abs(scale[dim - 1, dim - 1] / scale[j, j]) * ratio_delta
        out = max(out, scaled)
    return out


def companion_radius_for_nodes(nodes, derivs, delta: mp.mpf, dim: int, scaling: str) -> mp.mpf:
    scale = scaling_diag(nodes, scaling, dim)
    return max(companion_entry_radius(row, delta, dim, scale) for row in derivs)


def boundary_entry_radius(polys, active_modes, point, delta: mp.mpf, jet_order: int) -> mp.mpf:
    dim = jet_order - 1
    out = mp.mpf("0")
    for active_col in range(active_modes.cols):
        for state_col in range(dim):
            sens = mp.mpf("0")
            for k in range(1, jet_order):
                for jet_col in range(k):
                    n = k - 1 - jet_col
                    m = n - state_col
                    if m < 0 or m > n:
                        continue
                    coeff = mp.binomial(n, m) / mp.factorial(k)
                    projected = mp.fsum(
                        poly_derivative_value(poly, point, jet_col)
                        * active_modes[poly_idx, active_col]
                        for poly_idx, poly in enumerate(polys)
                    )
                    sens += abs(coeff * projected)
            out = max(out, delta * sens)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--eigenrow-json",
        default="endpoint_eigenrow_interval_propagation_200_consequence_theorem.json",
    )
    parser.add_argument(
        "--chart-ball-json",
        default="endpoint_grassmann_chart_ball_consequence_theorem.json",
    )
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--kind", choices=["kb"], default="kb")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--condition-matrix-order", type=int, default=200)
    parser.add_argument("--cache-dir", default=".endpoint_eigenrow_interval_cache")
    parser.add_argument("--active-reference-dps", type=int, default=90)
    parser.add_argument("--active-reference-matrix-order", type=int, default=90)
    parser.add_argument("--active-reference-cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
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
    parser.add_argument("--needed-trace-q", type=int, default=8)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
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
    parser.add_argument("--json-out", default="endpoint_coefficient_synchronized_200_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    eigenrow = load_json(args.eigenrow_json)
    chart_ball = load_json(args.chart_ball_json)
    derivative_radius = mp.mpf(str(eigenrow["maxDerivativeEntryRadius"]))
    synchronized = (
        int(eigenrow["conditionMatrixOrder"]) == int(eigenrow["controlledQuadratureOrder"])
        == args.condition_matrix_order
    )
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    dim = args.jet_order - 1
    max_q = args.needed_trace_q
    left_nodes = cheb_lobatto_nodes(args.order, left, s0)
    right_nodes = cheb_lobatto_nodes(args.order, s0, right)
    all_nodes = left_nodes + right_nodes

    print(
        "Endpoint coefficient synchronized certificate "
        f"order={args.order} matrix_order={args.condition_matrix_order}",
        flush=True,
    )
    cache, derivative_rows, min_gap, min_cos, hits = derivative_cache(args, all_nodes, max_q)
    left_derivs = [cache[node_key(node)] for node in left_nodes]
    right_derivs = [cache[node_key(node)] for node in right_nodes]
    e_s0 = nearest_derivs(cache, s0)
    e_left = nearest_derivs(cache, left)
    e_right = nearest_derivs(cache, right)

    reference_args = SimpleNamespace(**vars(args))
    reference_args.dps = args.active_reference_dps
    reference_args.matrix_order = args.active_reference_matrix_order
    reference_args.cache_dir = args.active_reference_cache_dir
    _ref_vals, e_reference_s0, _ref_lam = exact_derivative_at(
        reference_args,
        s0,
        max(args.max_trace_q, dim),
    )
    base, active_modes, active_idx = active_setup(args, e_reference_s0)
    c_left = active_endpoint_covectors(base["polys"], active_modes, left, e_left, args.jet_order)
    c_right = active_endpoint_covectors(base["polys"], active_modes, right, e_right, args.jet_order)

    left_companions = [
        scaling_diag(left_nodes, args.state_scaling, dim)
        * companion_matrix(row, dim)
        * (scaling_diag(left_nodes, args.state_scaling, dim) ** -1)
        for row in left_derivs
    ]
    right_companions = [
        scaling_diag(right_nodes, args.state_scaling, dim)
        * companion_matrix(row, dim)
        * (scaling_diag(right_nodes, args.state_scaling, dim) ** -1)
        for row in right_derivs
    ]
    companion_radius = max(
        companion_radius_for_nodes(left_nodes, left_derivs, derivative_radius, dim, args.state_scaling),
        companion_radius_for_nodes(right_nodes, right_derivs, derivative_radius, dim, args.state_scaling),
    )
    boundary_radius = max(
        boundary_entry_radius(base["polys"], active_modes, left, derivative_radius, args.jet_order),
        boundary_entry_radius(base["polys"], active_modes, right, derivative_radius, args.jet_order),
    )

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

    chart_cols = chart_ball["chartCoordinate"]
    threshold = mp.mpf(str(chart_ball["threshold"]))
    raw_capacity_row = raw_entry_radius_capacity(matrix_to_lists(endpoint), chart_cols, threshold)
    raw_capacity = mp.mpf(str(raw_capacity_row["uniformAbsoluteEntryRadius"]))
    left_cert = side_certificate(left_system, left_rhs, left_solution, left_scale, 0, dim, "left-to-s0 row flow")
    right_cert = side_certificate(right_system, right_rhs, right_solution, right_scale, args.order - 1, dim, "right-to-s0 row flow")
    coeff_cap, coeff_entry_radius, coeff_q = coeff_radius_capacity(left_cert, right_cert, mp.mpf("0"), raw_capacity)
    boundary_cap, boundary_entry_radius_at_cap, boundary_q = boundary_radius_capacity(left_cert, right_cert, mp.mpf("0"), raw_capacity)
    simultaneous_factor, simultaneous_radius, simultaneous_q = simultaneous_fraction_capacity(
        left_cert,
        right_cert,
        coeff_cap,
        boundary_cap,
        raw_capacity,
    )
    actual_endpoint_radius, actual_q = endpoint_entry_radius(
        left_cert,
        right_cert,
        companion_radius,
        boundary_radius,
    )
    actual_ratio_lower = raw_entry_ratio_lower(matrix_to_lists(endpoint), chart_cols, actual_endpoint_radius)
    input_pass = (
        companion_radius <= simultaneous_factor * coeff_cap
        and boundary_radius <= simultaneous_factor * boundary_cap
        and actual_q < 1
        and actual_endpoint_radius <= raw_capacity
    )

    data = {
        "theoremName": "endpoint coefficient synchronized 200 certificate",
        "order": args.order,
        "stateScaling": args.state_scaling,
        "conditionMatrixOrder": args.condition_matrix_order,
        "controlledQuadratureOrder": eigenrow["controlledQuadratureOrder"],
        "activeReferencePolicy": (
            "fixed legacy/base Krawczyk active source frame; synchronized "
            "200-point companion and boundary rows are projected into this "
            "same frame before applying the chart/Krawczyk budget"
        ),
        "activeReferenceDps": args.active_reference_dps,
        "activeReferenceMatrixOrder": args.active_reference_matrix_order,
        "centerSynchronized": bool(synchronized),
        "collocationNodeCount": len({node_key(node) for node in all_nodes}),
        "cacheHits": hits,
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "derivativeEntryRadius": f(derivative_radius),
        "derivativeEntryRadiusText": mp.nstr(derivative_radius, 16),
        "scaledCompanionAnalyticRadius": f(companion_radius),
        "scaledCompanionAnalyticRadiusText": mp.nstr(companion_radius, 16),
        "boundaryRowAnalyticRadius": f(boundary_radius),
        "boundaryRowAnalyticRadiusText": mp.nstr(boundary_radius, 16),
        "scaledCompanionMaxAbs": f(max(max_abs_matrix(left_companions), max_abs_matrix(right_companions))),
        "boundaryRowMaxAbs": f(max(matrix_max_abs(c_left), matrix_max_abs(c_right))),
        "chartCoordinate": chart_cols,
        "rawEndpointEntryCapacityForChartThreshold": f(raw_capacity),
        "rawEndpointEntryCapacityForChartThresholdText": mp.nstr(raw_capacity, 16),
        "rawEndpointRatioLowerAtActualRadius": f(actual_ratio_lower),
        "leftSystem": floatside(left_cert),
        "rightSystem": floatside(right_cert),
        "coefficientRadiusCapacityScaledCompanion": f(coeff_cap),
        "boundaryRadiusCapacityUnscaledRows": f(boundary_cap),
        "simultaneousCapacityFraction": f(simultaneous_factor),
        "simultaneousCoefficientRadius": f(simultaneous_factor * coeff_cap),
        "simultaneousBoundaryRadius": f(simultaneous_factor * boundary_cap),
        "actualEndpointEntryRadius": f(actual_endpoint_radius),
        "actualEndpointEntryRadiusText": mp.nstr(actual_endpoint_radius, 16),
        "actualKrawczykQ": f(actual_q),
        "scaledCompanionRadiusOverSimultaneousBudget": f(companion_radius / (simultaneous_factor * coeff_cap)),
        "boundaryRowRadiusOverSimultaneousBudget": f(boundary_radius / (simultaneous_factor * boundary_cap)),
        "endpointRadiusOverChartCapacity": f(actual_endpoint_radius / raw_capacity),
        "centerSynchronizationStatus": status(
            "200-point center synchronization",
            bool(synchronized),
            "Closed when the eigenrow interval and deterministic segment quadrature use the same matrix order.",
        ),
        "analyticCompanionInputStatus": status(
            "analytic scaled companion input below Krawczyk budget",
            bool(companion_radius <= simultaneous_factor * coeff_cap),
            "Closed when the companion-map perturbation radius fits inside the simultaneous Krawczyk coefficient budget.",
        ),
        "analyticBoundaryInputStatus": status(
            "analytic boundary-row input below Krawczyk budget",
            bool(boundary_radius <= simultaneous_factor * boundary_cap),
            "Closed when the Green-concomitant boundary-row perturbation radius fits inside the simultaneous Krawczyk boundary budget.",
        ),
        "synchronizedKrawczykInputStatus": status(
            "synchronized endpoint Krawczyk input",
            bool(synchronized and input_pass),
            "Closed when the synchronized center is used and both analytic input radii keep the endpoint map inside the persistent Pluecker chart.",
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint coefficient synchronized 200 certificate")
    print(f"  derivative radius: {mp.nstr(derivative_radius, 12)}")
    print(f"  companion radius: {mp.nstr(companion_radius, 12)} / simultaneous {mp.nstr(simultaneous_factor * coeff_cap, 12)}")
    print(f"  boundary radius: {mp.nstr(boundary_radius, 12)} / simultaneous {mp.nstr(simultaneous_factor * boundary_cap, 12)}")
    print(f"  endpoint radius: {mp.nstr(actual_endpoint_radius, 12)} / chart cap {mp.nstr(raw_capacity, 12)}")
    print(f"  actual Krawczyk q: {mp.nstr(actual_q, 12)}")
    print(f"  synchronized input closed: {synchronized and input_pass}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
