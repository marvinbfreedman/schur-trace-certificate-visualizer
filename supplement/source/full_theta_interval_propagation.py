#!/usr/bin/env python3
r"""Interval propagation of the omitted full-theta tail through Psi -> E_u.

This is the final formal step after the finite source-tail certificate.  The
finite script proved that

    S_{Phi<=8} - S_{tilde Phi_3}

is tiny in the source/Riesz coordinates.  Here we propagate the omitted theta
profile

    R(v) = sum_{n>=9} phi_n(v)

through the normalized map

    Psi -> h_u^Psi(s) = K_red^Psi(s,u)
         -> E_u^Psi f = (B_P[h_u^Psi,f](s0), (P^*h_u^Psi)(s0) f(s0)).

The propagation is intentionally conservative.  For each derivative order
0<=k<=8, it encloses R^(k)(v) on the whole source window by an interval
[-tau_k,tau_k].  Those intervals are inserted into the Taylor coefficients of
Psi in every normalized ratio

    Psi(x+r)/Psi(x).

The resulting interval rows enclose Delta E = E_{Phi}-E_{Phi<=8} on the finite
source quadrature grid.  The source Gram perturbation is bounded by

    ||Delta S|| <= 2 ||E_8|| ||Delta E|| + ||Delta E||^2.

The output is a finite-rank interval certificate in the same A-normalized
source/Riesz coordinates as ``source_side_riesz_rank_theorem.py``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders, submatrix_rows, sym  # noqa: E402
from full_theta_source_tail_certificate import (  # noqa: E402
    finite_riesz_margin,
    full_weights,
    h_derivatives_for_pieces,
    pieces_from_mode_weights,
    psi_series,
    weighted_source_matrix_for_pieces,
)
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, source_nodes  # noqa: E402
from lagrange_energy_control_certificate import boundary_functional_row, eval_functional_row  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, endpoint_b_quadrature  # noqa: E402
from source_quadrature_interval_envelope import (  # noqa: E402
    adjoint_source_value_interval,
    boundary_functional_row_interval,
    div_any,
    eval_functional_row_interval,
    interval_abs_sup,
    mul_any,
    sum_any,
    trace_green_concomitant_row_interval,
)
from source_side_quadrature_refinement import weight_source_rows, weights_for_nodes  # noqa: E402
from source_side_riesz_rank_theorem import op_norm_rect  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value, trace_green_concomitant_row  # noqa: E402


def iv_ball(radius):
    radius = mp.mpf(radius)
    return mp.iv.mpf([-radius, radius])


def iv_point(x):
    x = mp.mpf(x)
    return mp.iv.mpf([x, x])


def scale_any(a, c):
    return [x * c for x in a]


def mode_derivative_bound(n, k, vmin):
    n2 = n * n
    c = mp.pi * n2
    y = mp.e**vmin
    total = mp.mpf("0")
    for coeff, beta in ((4 * c * c, mp.mpf("2.25")), (-6 * c, mp.mpf("1.25"))):
        base = abs(coeff) * y**beta * mp.e ** (-c * y)
        factor = 1 + beta + c * y
        total += base * factor**k
    return total


def infinite_tail_derivative_bounds(start_mode, vmin, max_order):
    bounds = []
    ratios = []
    for k in range(max_order + 1):
        first = mode_derivative_bound(start_mode, k, vmin)
        second = mode_derivative_bound(start_mode + 1, k, vmin)
        q = second / first if first else mp.mpf("0")
        if q >= 1:
            raise ValueError(f"tail ratio is not contracting for k={k}: {q}")
        bounds.append(first / (1 - q))
        ratios.append(q)
    return bounds, ratios


def psi_series_with_tail_interval(v0, pieces, tail_bounds, n):
    finite = psi_series(v0, pieces, n)
    out = []
    for k in range(n):
        out.append(iv_point(finite[k]) + iv_ball(tail_bounds[k] / mp.factorial(k)))
    return out


def ratio_series_with_tail_interval(x0, r, pieces, tail_bounds, n):
    numerator = psi_series_with_tail_interval(x0 + r, pieces, tail_bounds, n)
    denominator = psi_series_with_tail_interval(x0, pieces, tail_bounds, n)
    return div_any(numerator, denominator), denominator[0]


def ratio_value_with_tail_interval(x0, r, pieces, tail_bounds):
    numerator = psi_series_with_tail_interval(x0 + r, pieces, tail_bounds, 1)[0]
    denominator = psi_series_with_tail_interval(x0, pieces, tail_bounds, 1)[0]
    return numerator / denominator, denominator


def ratio_series_finite(x0, r, pieces, n):
    return div_any(psi_series(x0 + r, pieces, n), psi_series(x0, pieces, n))


def ratio_value_finite(x0, r, pieces):
    return psi_series(x0 + r, pieces, 1)[0] / psi_series(x0, pieces, 1)[0]


def weight_series(center, omega, n):
    # Series in delta_s for (r+(s+u)/2) cosh(omega(r+(s+u)/2)).
    y = [mp.mpf("0") for _ in range(n)]
    y[0] = center
    if n > 1:
        y[1] = mp.mpf("0.5")
    ep = exp_series(scale_any(y, omega))
    em = exp_series(scale_any(y, -omega))
    cosh = [(ep[k] + em[k]) * mp.mpf("0.5") for k in range(n)]
    return mul_any(y, cosh)


def exp_series(a):
    n = len(a)
    out = [mp.mpf("0") for _ in range(n)]
    out[0] = mp.e**a[0]
    for i in range(1, n):
        out[i] = mp.fsum(k * a[k] * out[i - k] for k in range(1, i + 1)) / i
    return out


def build_s_ratio_cache(args, pieces8, tail_bounds, r_nodes):
    n = args.jet_order
    s0 = mp.mpf(args.s0)
    cache = []
    for r in r_nodes:
        ratio_iv, denom_s = ratio_series_with_tail_interval(s0, r, pieces8, tail_bounds, n)
        ratio_fin = ratio_series_finite(s0, r, pieces8, n)
        denom_lower = min(abs(mp.mpf(denom_s.a)), abs(mp.mpf(denom_s.b)))
        cache.append((ratio_iv, ratio_fin, denom_lower))
    return cache


def h_delta_derivatives_interval(args, pieces8, tail_bounds, u, r_nodes, r_weights, s_cache):
    n = args.jet_order
    s0 = mp.mpf(args.s0)
    omega = mp.mpf(args.omega)
    coeffs = [mp.iv.mpf([0, 0]) for _ in range(n)]
    finite_coeffs = [mp.mpf("0") for _ in range(n)]
    denom_abs_lowers = []

    for r, weight, cached in zip(r_nodes, r_weights, s_cache):
        ratio_s, ratio_s_fin, denom_s_lower = cached
        ratio_u, denom_u = ratio_value_with_tail_interval(u, r, pieces8, tail_bounds)
        ratio_u_fin = ratio_value_finite(u, r, pieces8)
        denom_abs_lowers.extend(
            [
                denom_s_lower,
                min(abs(mp.mpf(denom_u.a)), abs(mp.mpf(denom_u.b))),
            ]
        )
        w = weight_series(r + (s0 + u) / 2, omega, n)
        term = scale_any(mul_any(w, ratio_s), ratio_u * weight)
        finite_term = scale_any(mul_any(w, ratio_s_fin), ratio_u_fin * weight)
        for k in range(n):
            coeffs[k] += term[k]
            finite_coeffs[k] += finite_term[k]

    derivs = []
    for k in range(n):
        derivs.append((coeffs[k] - finite_coeffs[k]) * mp.factorial(k))
    return derivs, min(denom_abs_lowers)


def source_delta_rows_interval(args, polys, e_derivs, pieces8, tail_bounds, r_nodes, r_weights, s_cache, u):
    h_delta, denom_lower = h_delta_derivatives_interval(
        args,
        pieces8,
        tail_bounds,
        u,
        r_nodes,
        r_weights,
        s_cache,
    )
    pstar = adjoint_source_value_interval(e_derivs, h_delta, args.jet_order)
    brow = trace_green_concomitant_row_interval(e_derivs, h_delta, args.jet_order)
    brow_poly = boundary_functional_row_interval(polys, mp.mpf(args.s0), brow)
    eval_poly = eval_functional_row_interval(polys, mp.mpf(args.s0), pstar)
    return [brow_poly, eval_poly], h_delta, denom_lower


def interval_scaled_rows(args, base, e_derivs, pieces8, tail_bounds):
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.source_rmax), args.source_order)
    s_cache = build_s_ratio_cache(args, pieces8, tail_bounds, r_nodes)
    nodes = source_nodes(args)
    weights = weights_for_nodes(nodes)
    rows = []
    max_h_bounds = [mp.mpf("0") for _ in range(args.jet_order)]
    min_denom_lower = mp.inf
    row_norm_bounds = []

    for block, u in enumerate(nodes):
        source_rows, h_delta, denom_lower = source_delta_rows_interval(
            args,
            base["polys"],
            e_derivs,
            pieces8,
            tail_bounds,
            r_nodes,
            r_weights,
            s_cache,
            u,
        )
        min_denom_lower = min(min_denom_lower, denom_lower)
        for k, value in enumerate(h_delta):
            max_h_bounds[k] = max(max_h_bounds[k], interval_abs_sup(value))
        scale = mp.sqrt(weights[block])
        for local_row in source_rows:
            scaled = []
            for j in range(base["scaledModes"].cols):
                value = sum_any(
                    local_row[k] * base["scaledModes"][k, j]
                    for k in range(len(local_row))
                )
                scaled.append(value * scale)
            rows.append(scaled)
            row_norm_bounds.append(mp.sqrt(mp.fsum(interval_abs_sup(v) ** 2 for v in scaled)))
        if (block + 1) % max(1, len(nodes) // 5) == 0:
            print(f"  interval rows {block + 1}/{len(nodes)}", flush=True)
    return nodes, weights, rows, max_h_bounds, min_denom_lower, row_norm_bounds


def interval_matrix_frobenius_bound(rows):
    total = mp.mpf("0")
    for row in rows:
        for value in row:
            total += interval_abs_sup(value) ** 2
    return mp.sqrt(total)


def finite_response_margin(args, base, source_gram):
    vals, vecs = mp.eigsy(sym(source_gram), eigvals_only=False)
    active_idx = [len(vals) - 2, len(vals) - 1]
    active = columns(vecs, active_idx)
    gap = vals[-2] - vals[-3]
    _lmat, finite_svals, rank_margin, response_norm = finite_riesz_margin(args, base, active)
    return vals, active_idx, gap, finite_svals, rank_margin, response_norm


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=65)
    parser.add_argument("--full-nmax", type=int, default=8)
    parser.add_argument("--tail-start", type=int, default=9)
    parser.add_argument("--source-order", type=int, default=16)
    parser.add_argument("--source-rmax", default="10")
    parser.add_argument("--max-deriv", type=int, default=8)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-6")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--iv-dps", type=int, default=40)
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
    parser.add_argument("--json-out", default="full_theta_interval_propagation.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps

    vmin = min(mp.mpf(args.source_min), mp.mpf(args.s0))
    tail_bounds, tail_ratios = infinite_tail_derivative_bounds(
        args.tail_start,
        vmin,
        args.jet_order - 1,
    )
    pieces8 = pieces_from_mode_weights(full_weights(args.full_nmax))

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace/high block basis={args.basis} max_q={max_q}", flush=True)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)
    base_args = SimpleNamespace(**vars(args))
    base_args.source_grid = 9
    base = local_case(base_args, args.basis, e_derivs_at_s0)

    print("Building finite E_8 source Gram", flush=True)
    _nodes, _weights, e8, s8 = weighted_source_matrix_for_pieces(
        args,
        base,
        e_derivs_at_s0,
        pieces8,
    )
    e8_norm = op_norm_rect(e8)
    s8_norm = op_norm_rect(e8) ** 2
    svals, active_idx, gap, finite_svals, rank_margin, response_norm = finite_response_margin(args, base, s8)

    print("Propagating omitted tail intervals through E_u", flush=True)
    nodes, weights, delta_rows, max_h_bounds, min_denom_lower, row_norm_bounds = interval_scaled_rows(
        args,
        base,
        e_derivs_at_s0,
        pieces8,
        tail_bounds,
    )
    delta_e_bound = interval_matrix_frobenius_bound(delta_rows)
    delta_s_bound = 2 * e8_norm * delta_e_bound + delta_e_bound**2
    alpha = 4 * delta_s_bound / gap if gap else mp.inf
    lower_after_tail = rank_margin * (1 - alpha) - response_norm * alpha
    pass_gap = bool(delta_s_bound < gap / 4)
    pass_rank = bool(pass_gap and lower_after_tail > 0)

    print(f"  max ||Delta E row||_F bound={fmt(delta_e_bound, 12)}", flush=True)
    print(f"  ||Delta S|| bound={fmt(delta_s_bound, 12)}", flush=True)
    print(f"  ||Delta S||/||S8||={fmt(delta_s_bound/s8_norm if s8_norm else mp.inf, 12)}", flush=True)
    print(f"  gap={fmt(gap, 12)} alpha={fmt(alpha, 12)}", flush=True)
    print(f"  lower_after_tail={fmt(lower_after_tail, 12)} pass={pass_rank}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "sourceGrid": args.source_grid,
        "fullNmax": args.full_nmax,
        "tailStart": args.tail_start,
        "sourceOrder": args.source_order,
        "sourceRmax": f(mp.mpf(args.source_rmax)),
        "vMin": f(vmin),
        "ivDps": args.iv_dps,
        "tailDerivativeBounds": [f(x) for x in tail_bounds],
        "tailGeometricRatios": [f(x) for x in tail_ratios],
        "maxPropagatedHDerivativeBounds": [f(x) for x in max_h_bounds],
        "minDenominatorAbsLower": f(min_denom_lower),
        "sourceNodesFirstLast": [f(nodes[0]), f(nodes[-1])],
        "sourceWeightsFirstLast": [f(weights[0]), f(weights[-1])],
        "rowNormBounds": [f(x) for x in row_norm_bounds],
        "deltaEOperatorBound": f(delta_e_bound),
        "finiteE8OperatorNorm": f(e8_norm),
        "finiteS8OperatorNorm": f(s8_norm),
        "deltaSOperatorBound": f(delta_s_bound),
        "deltaSRelativeToS8": f(delta_s_bound / s8_norm if s8_norm else mp.inf),
        "activeEigenvaluesS8": [f(svals[i]) for i in active_idx],
        "complementTopEigenvalueS8": f(svals[-3]),
        "spectralGapToComplementS8": f(gap),
        "tailProjectorAlpha": f(alpha),
        "tailGapConditionPasses": pass_gap,
        "finiteActiveResponseSingularValuesS8": [f(x) for x in finite_svals],
        "finiteRankMarginS8": f(rank_margin),
        "responseOperatorNorm": f(response_norm),
        "lowerBoundAfterOmittedTail": f(lower_after_tail),
        "rankPersistsForFullPhi": pass_rank,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Interval propagation of the omitted n>=tailStart theta profile through "
            "the normalized source-row map Psi -> E_u^Psi on the finite source grid. "
            "The bound uses derivative envelopes for the theta tail, interval quotient "
            "arithmetic for Psi(x+r)/Psi(x), and ||Delta S|| <= 2||E8||||Delta E||+||Delta E||^2."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
