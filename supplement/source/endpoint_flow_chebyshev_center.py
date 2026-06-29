#!/usr/bin/env python3
r"""Chebyshev-collocation center for the endpoint adjoint flow.

Target system:

    Y'(s) = A(s)Y(s),        P^*K=0,

where ``A(s)`` is the first-order companion matrix of the formal adjoint trace
operator.  Unlike the older endpoint scripts, this center construction uses:

* exact local confluent eigen-derivatives for the trace row ``e(s)`` at
  Chebyshev nodes;
* global Chebyshev-Lobatto collocation on ``[left,s0]`` and ``[s0,right]``;
* no finite differences of the moving eigenrow and no piecewise-constant
  ``expm`` propagation.

The output is still a numerical center, not an interval proof.  Its purpose is
to test whether endpoint maps stabilize under spectral refinement before
building a validated/ball-arithmetic version.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    companion_matrix,
    endpoint_active_row,
    endpoint_map,
    left_obstruction,
    mat_frobenius,
    matrix_rank_gram,
    min_norm_solve,
    row_norm,
    status,
)
from adjoint_green_jump_conditions import jump_matrix, solve_jump  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def vec_dot(a, b) -> mp.mpf:
    return mp.fsum(a[i] * b[i] for i in range(len(a)))


def flip_derivs(e_derivs):
    return [[-value for value in row] for row in e_derivs]


def cheb_lobatto_nodes(order: int, left: mp.mpf, right: mp.mpf) -> list[mp.mpf]:
    if order < 2:
        raise SystemExit("Chebyshev order must be at least 2")
    mid = (left + right) / 2
    half = (right - left) / 2
    n = order - 1
    return [mid + half * mp.cos(mp.pi * j / n) for j in range(order)]


def cheb_lobatto_diff_matrix(order: int, left: mp.mpf, right: mp.mpf) -> mp.matrix:
    """Return differentiation matrix d/ds at Lobatto nodes cos(pi*j/N)."""
    n = order - 1
    xs = [mp.cos(mp.pi * j / n) for j in range(order)]
    c = [mp.mpf("2") if j in (0, n) else mp.mpf("1") for j in range(order)]
    d = mp.matrix(order)
    for i in range(order):
        for j in range(order):
            if i == j:
                continue
            d[i, j] = (c[i] / c[j]) * ((-1) ** (i + j)) / (xs[i] - xs[j])
    for i in range(order):
        d[i, i] = -mp.fsum(d[i, j] for j in range(order) if j != i)
    return (mp.mpf("2") / (right - left)) * d


def node_key(value: mp.mpf) -> str:
    return mp.nstr(value, 80)


def exact_derivative_at(args, s: mp.mpf, max_q: int):
    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    cache_key = {
        "s": mp.nstr(s, 80),
        "max_q": max_q,
        "jet_order": args.jet_order,
        "dps": args.dps,
        "matrix_order": args.matrix_order,
        "matrix_rmax": args.matrix_rmax,
    }
    cache_path = None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256(
            json.dumps(cache_key, sort_keys=True).encode("utf-8")
        ).hexdigest()[:24]
        cache_path = cache_dir / f"trace_derivs_{digest}.json"
        if cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            vals = [mp.mpf(value) for value in cached["vals"]]
            e_derivs = [
                [mp.mpf(value) for value in row]
                for row in cached["eDerivs"]
            ]
            lam_derivs = [mp.mpf(value) for value in cached["lambdaDerivs"]]
            return vals, e_derivs, lam_derivs

    local_args = SimpleNamespace(**vars(args))
    local_args.s0 = str(s)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(local_args, max_q)
    if cache_path is not None:
        cache_path.write_text(
            json.dumps(
                {
                    "key": cache_key,
                    "vals": [mp.nstr(value, 90) for value in vals],
                    "eDerivs": [
                        [mp.nstr(value, 90) for value in row]
                        for row in e_derivs
                    ],
                    "lambdaDerivs": [mp.nstr(value, 90) for value in lam_derivs],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return vals, e_derivs, lam_derivs


def build_derivative_cache(args, nodes: list[mp.mpf], max_q: int, full_q: int):
    unique = sorted({node_key(node): node for node in nodes}.values())
    cache = {}
    previous = None
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    rows = []
    for idx, node in enumerate(unique):
        vals, e_derivs, lam_derivs = exact_derivative_at(args, node, max_q)
        if previous is not None:
            cos = vec_dot(previous, e_derivs[0])
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
                "lambda0": f(vals[0]),
                "lambda1": f(vals[1]),
                "gap": f(gap),
                "e8": f(e_derivs[0][-1]),
            }
        )

    s0 = mp.mpf(args.s0)
    vals_full, e_full, lam_full = exact_derivative_at(args, s0, full_q)
    s0_short = cache.get(node_key(s0))
    if s0_short is not None and vec_dot(s0_short[0], e_full[0]) < 0:
        e_full = flip_derivs(e_full)
    return cache, e_full, rows, min_gap, min_cos


def scaling_diag(nodes: list[mp.mpf], mode: str, dim: int) -> mp.matrix:
    scale = mp.matrix(dim, dim)
    if mode == "raw":
        for j in range(dim):
            scale[j, j] = 1
        return scale
    if mode != "taylor":
        raise ValueError(mode)
    radius = (max(nodes) - min(nodes)) / 2
    for j in range(dim):
        scale[j, j] = radius**j / mp.factorial(j)
    return scale


def scaled_companion(amat: mp.matrix, scale: mp.matrix) -> mp.matrix:
    return scale * amat * (scale ** -1)


def flow_collocation(
    nodes: list[mp.mpf],
    derivs,
    init_index: int,
    target_index: int,
    scaling: str,
) -> mp.matrix:
    order = len(nodes)
    dim = 8
    left = min(nodes)
    right = max(nodes)
    dmat = cheb_lobatto_diff_matrix(order, left, right)
    scale = scaling_diag(nodes, scaling, dim)
    scale_inv = scale ** -1
    amats = [
        scaled_companion(companion_matrix(derivs[j], dim), scale)
        for j in range(order)
    ]
    size = order * dim
    system = mp.matrix(size, size)
    rhs = mp.matrix(size, dim)

    def idx(node: int, comp: int) -> int:
        return node * dim + comp

    for j in range(order):
        for m in range(dim):
            row = idx(j, m)
            if j == init_index:
                system[row, idx(j, m)] = 1
                rhs[row, m] = 1
                continue
            for k in range(order):
                system[row, idx(k, m)] += dmat[j, k]
            amat = amats[j]
            for ell in range(dim):
                system[row, idx(j, ell)] -= amat[m, ell]

    sol = mp.matrix(size, dim)
    for col in range(dim):
        rhs_col = mp.matrix(size, 1)
        for row in range(size):
            rhs_col[row] = rhs[row, col]
        sol_col = mp.lu_solve(system, rhs_col)
        for row in range(size):
            sol[row, col] = sol_col[row]
    out = mp.matrix(dim, dim)
    for m in range(dim):
        for col in range(dim):
            out[m, col] = sol[idx(target_index, m), col]
    return scale_inv * out * scale


def source_endpoint_rows(args, base, active_modes, e_left, e_right, flow_right, e_s0):
    s0 = mp.mpf(args.s0)
    left = mp.mpf(args.constraint_min)
    right = mp.mpf(args.constraint_max)
    jmat = jump_matrix(e_s0, args.jet_order - 1)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    zero_left = mp.matrix(args.jet_order - 1, 1)
    rows = []
    for t in t_values:
        h_derivs = source_derivatives(mp.pi, s0, t, args.jet_order, r_nodes, r_weights)
        boundary_source = trace_green_concomitant_row(e_s0, h_derivs, args.jet_order)
        pstar = adjoint_source_value(e_s0, h_derivs, args.jet_order)
        eval_source = [pstar] + [mp.mpf("0") for _ in range(args.jet_order - 2)]
        for component, source_row in (("boundary", boundary_source), ("adjointEval", eval_source)):
            delta, jump_residual, _target = solve_jump(jmat, source_row)
            delta_vec = mp.matrix(len(delta), 1)
            for i, value in enumerate(delta):
                delta_vec[i] = value
            b_values = endpoint_active_row(
                base["polys"],
                active_modes,
                e_left,
                e_right,
                zero_left,
                flow_right * delta_vec,
                left,
                right,
            )
            rows.append(
                {
                    "sourceNode": f(t),
                    "component": component,
                    "jumpNorm": f(row_norm(delta)),
                    "endpointVector": [f(value) for value in b_values],
                    "endpointBeforeNorm": f(row_norm(b_values)),
                    "jumpResidualNorm": f(row_norm(jump_residual)),
                }
            )
    return rows


def build_center(args, order: int) -> dict:
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    left_nodes = cheb_lobatto_nodes(order, left, s0)
    right_nodes = cheb_lobatto_nodes(order, s0, right)
    all_nodes = left_nodes + right_nodes
    cache, e_s0, derivative_rows, min_gap, min_cos = build_derivative_cache(
        args,
        all_nodes,
        args.jet_order - 1,
        max(args.max_trace_q, args.jet_order - 1),
    )

    left_derivs = [cache[node_key(node)] for node in left_nodes]
    right_derivs = [cache[node_key(node)] for node in right_nodes]
    # Lobatto order is high->low.  On [left,s0], init at s0 is index 0 and
    # target left is index order-1.  On [s0,right], init at s0 is index order-1
    # and target right is index 0.
    flow_left = flow_collocation(
        left_nodes,
        left_derivs,
        0,
        order - 1,
        args.state_scaling,
    )
    flow_right = flow_collocation(
        right_nodes,
        right_derivs,
        order - 1,
        0,
        args.state_scaling,
    )

    base, active_modes, active_idx = active_setup(args, e_s0)
    e_left = left_derivs[-1]
    e_right = right_derivs[0]
    endpoint_mat, _zero = endpoint_map(
        base["polys"],
        active_modes,
        e_left,
        e_right,
        flow_left,
        flow_right,
        left,
        right,
    )
    rank, eigs, rank_threshold = matrix_rank_gram(
        endpoint_mat * endpoint_mat.T,
        mp.mpf(args.trace_tol),
    )
    obstruction, obstruction_eigs = left_obstruction(endpoint_mat)
    obstruction_residual = endpoint_mat.T * obstruction
    obstruction_residual_norm = row_norm(
        [obstruction_residual[i] for i in range(obstruction_residual.rows)]
    )
    raw_rows = source_endpoint_rows(
        args,
        base,
        active_modes,
        e_left,
        e_right,
        flow_right,
        e_s0,
    )
    max_residual = mp.mpf("0")
    max_z_norm = mp.mpf("0")
    max_obstruction_pairing = mp.mpf("0")
    max_obstruction_relative_pairing = mp.mpf("0")
    rows = []
    for row in raw_rows:
        b_values = [mp.mpf(str(value)) for value in row["endpointVector"]]
        z, endpoint_residual, _gram = min_norm_solve(
            endpoint_mat,
            [-value for value in b_values],
        )
        residual_norm = row_norm([endpoint_residual[i] for i in range(endpoint_residual.rows)])
        z_norm = row_norm([z[i] for i in range(z.rows)])
        beta_before = row_norm(b_values)
        obstruction_pairing = mp.fsum(
            obstruction[i] * b_values[i] for i in range(len(b_values))
        )
        obstruction_relative = abs(obstruction_pairing) / max(mp.mpf("1"), beta_before)
        max_residual = max(max_residual, residual_norm)
        max_z_norm = max(max_z_norm, z_norm)
        max_obstruction_pairing = max(max_obstruction_pairing, abs(obstruction_pairing))
        max_obstruction_relative_pairing = max(
            max_obstruction_relative_pairing,
            obstruction_relative,
        )
        row.update(
            {
                "leftObstructionPairing": f(obstruction_pairing),
                "leftObstructionRelativePairing": f(obstruction_relative),
                "selectedZNorm": f(z_norm),
                "endpointResidualNorm": f(residual_norm),
                "endpointRelativeResidual": f(residual_norm / max(mp.mpf("1"), beta_before)),
                "selectedZ": [f(z[i]) for i in range(z.rows)],
            }
        )
        rows.append(row)

    sigmas = [mp.sqrt(max(mp.mpf("0"), value)) for value in eigs]
    return {
        "chebyshevOrder": order,
        "stateScaling": args.state_scaling,
        "activeDim": len(active_idx),
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min_gap),
        "traceFieldMinConsecutiveCos": f(min_cos),
        "endpointMapRank": rank,
        "endpointMapGramEigenvalues": [f(value) for value in eigs],
        "endpointMapRankThreshold": f(rank_threshold),
        "sigmaMin": f(min(sigmas)),
        "sigmaMax": f(max(sigmas)),
        "endpointMap": [
            [f(endpoint_mat[i, j]) for j in range(endpoint_mat.cols)]
            for i in range(endpoint_mat.rows)
        ],
        "leftObstructionVector": [f(obstruction[i]) for i in range(obstruction.rows)],
        "leftObstructionGramEigenvalues": [f(value) for value in obstruction_eigs],
        "leftObstructionResidualNorm": f(obstruction_residual_norm),
        "maxLeftObstructionPairing": f(max_obstruction_pairing),
        "maxLeftObstructionRelativePairing": f(max_obstruction_relative_pairing),
        "endpointMapFrobenius": f(mat_frobenius(endpoint_mat)),
        "flowLeftFrobenius": f(mat_frobenius(flow_left)),
        "flowRightFrobenius": f(mat_frobenius(flow_right)),
        "maxEndpointResidualNorm": f(max_residual),
        "maxSelectedZNorm": f(max_z_norm),
        "derivativeRows": derivative_rows,
        "rows": rows,
        "statuses": {
            "activeEndpointMapFullRowRank": status(
                "active endpoint map full row rank",
                rank == len(active_idx),
                "Chebyshev-collocation center has full sampled active row rank.",
            ),
            "actualEndpointRhsInRange": status(
                "actual endpoint RHS in Range(M)",
                max_residual < mp.mpf("1e-20"),
                "For the actual source rows, the endpoint vectors b(d) are in the computed endpoint-map range.",
            ),
            "activeEndpointKilled": status(
                "active endpoint concomitant killed",
                max_residual < mp.mpf("1e-20"),
                "Minimum-norm z solves Mz=-b(d) on the active two-dimensional source space.",
            ),
            "finiteLeftObstructionAnnihilatesActualRows": status(
                "finite left obstruction annihilates actual rows",
                max_obstruction_relative_pairing < mp.mpf("1e-20"),
                "The resolved left endpoint obstruction annihilates the actual source endpoint vectors.",
            ),
        },
    }


def matrix_from_rows(rows: list[list[float]]) -> mp.matrix:
    return mp.matrix([[mp.mpf(str(value)) for value in row] for row in rows])


def compare_centers(left: dict, right: dict, margin: mp.mpf) -> dict:
    a = matrix_from_rows(left["endpointMap"])
    b = matrix_from_rows(right["endpointMap"])
    diff = b - a
    diff_norm = mat_frobenius(diff)
    return {
        "leftOrder": left["chebyshevOrder"],
        "rightOrder": right["chebyshevOrder"],
        "frobeniusDifference": f(diff_norm),
        "rankMargin": f(margin),
        "differenceOverMargin": f(diff_norm / margin),
        "passesRankMargin": bool(diff_norm < margin),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", default="5 7")
    parser.add_argument("--rank-margin", default="2.77565955249e6")
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
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
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
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
    parser.add_argument("--json-out", default="endpoint_flow_chebyshev_center.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    orders = parse_ints(args.orders)
    if len(orders) < 1:
        raise SystemExit("--orders must contain at least one order")
    rank_margin = mp.mpf(args.rank_margin)

    centers = []
    for order in orders:
        print(f"Building Chebyshev endpoint center order={order}", flush=True)
        center = build_center(args, order)
        centers.append(center)
        print(
            f"  order={order} rank={center['endpointMapRank']} "
            f"||M||_F={fmt(mp.mpf(center['endpointMapFrobenius']), 12)} "
            f"sigma_min={fmt(mp.mpf(center['sigmaMin']), 12)} "
            f"endpoint_res={fmt(mp.mpf(center['maxEndpointResidualNorm']), 8)}",
            flush=True,
        )

    comparisons = [
        compare_centers(left, right, rank_margin)
        for left, right in zip(centers[:-1], centers[1:])
    ]
    stable = bool(comparisons and all(row["passesRankMargin"] for row in comparisons))
    data = {
        "theoremName": "endpoint flow Chebyshev center",
        "interval": [f(mp.mpf(args.constraint_min)), f(mp.mpf(args.constraint_max))],
        "s0": f(mp.mpf(args.s0)),
        "orders": orders,
        "stateScaling": args.state_scaling,
        "rankMargin": f(rank_margin),
        "centers": centers,
        "comparisons": comparisons,
        "chebyshevCenterStabilityStatus": status(
            "Chebyshev endpoint center stability",
            stable,
            (
                "Closed only if successive Chebyshev endpoint maps differ by "
                "less than the rank margin.  This is still a numerical center "
                "test, not an interval proof."
            ),
        ),
        "nextTarget": (
            "If the center stabilizes, wrap the Chebyshev collocation solve in "
            "ball arithmetic/Taylor remainder bounds to prove the endpoint-map "
            "rank margin.  If it does not stabilize, return to analytic "
            "obstruction annihilation."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint flow Chebyshev center")
    for row in comparisons:
        print(
            f"  {row['leftOrder']}->{row['rightOrder']} "
            f"diff={row['frobeniusDifference']:.12e} "
            f"diff/margin={row['differenceOverMargin']:.12e} "
            f"pass={row['passesRankMargin']}",
            flush=True,
        )
    print(f"  stable below rank margin: {stable}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
