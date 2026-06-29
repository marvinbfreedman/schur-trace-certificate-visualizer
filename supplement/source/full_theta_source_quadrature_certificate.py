#!/usr/bin/env python3
r"""Continuum source-quadrature certificate for the full-theta source rows.

This promotes the source-side noncollapse certificate from a weighted source
grid to the continuum source window.

For the finite full theta profile Psi_8 = Phi_{<=8}, we build

    h_u(s) = K_red^{Psi_8}(s,u),
    E_u f = (B_P[h_u,f](s0), (P^*h_u)(s0) f(s0)).

Then

    S_8 = int_a^b E_u^* E_u du

is compared with the composite trapezoid source matrix S_{8,h}.  The
trapezoid error is bounded by

    ||S_8-S_{8,h}|| <= (b-a) h^2/12 sup_u ||d_u^2(E_u^*E_u)||.

The supremum is enclosed by real interval arithmetic over a partition of the
source window.  Finally the omitted n>=9 interval-propagation certificate is
added to the quadrature error before applying the Riesz projector margin.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders, sym  # noqa: E402
from full_theta_source_tail_certificate import (  # noqa: E402
    finite_riesz_margin,
    full_weights,
    pieces_from_mode_weights,
    weighted_source_matrix_for_pieces,
)
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from lagrange_energy_control_certificate import boundary_functional_row, eval_functional_row  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, endpoint_b_quadrature  # noqa: E402
from source_quadrature_interval_envelope import (  # noqa: E402
    adjoint_source_value_interval,
    boundary_functional_row_interval,
    div_any,
    eval_functional_row_interval,
    interval_abs_sup,
    interval_matrix_norm_bound,
    mul_any,
    sum_any,
    trace_green_concomitant_row_interval,
)
from source_side_riesz_rank_theorem import op_norm_rect  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value, trace_green_concomitant_row  # noqa: E402


def zero(n):
    return [mp.mpf("0") for _ in range(n)]


def is_interval(x):
    return hasattr(x, "_mpi_")


def scalar_add(x, y):
    if is_interval(y) and not is_interval(x):
        return y + x
    return x + y


def scalar_mul(x, y):
    if is_interval(y) and not is_interval(x):
        return y * x
    return x * y


def add(a, b):
    return [scalar_add(x, y) for x, y in zip(a, b)]


def scale(a, k):
    return [scalar_mul(x, k) for x in a]


def exp_any(x):
    return mp.iv.exp(x) if hasattr(x, "_mpi_") else mp.e**x


def exp_series_any(a):
    n = len(a)
    out = zero(n)
    out[0] = exp_any(a[0])
    for i in range(1, n):
        out[i] = sum_any(
            scalar_mul(scalar_mul(a[k], k), out[i - k]) for k in range(1, i + 1)
        ) / i
    return out


def psi_series_any(v0, pieces, n):
    ev0 = exp_any(v0)
    e_delta = [mp.mpf(1) / mp.factorial(k) for k in range(n)]
    out = zero(n)
    for coeff, beta, c in pieces:
        exponent = scale(e_delta, scalar_mul(ev0, -c))
        exponent[0] = scalar_add(exponent[0], scalar_mul(v0, beta))
        if n > 1:
            exponent[1] = scalar_add(exponent[1], beta)
        out = add(out, scale(exp_series_any(exponent), coeff))
    return out


def ratio_series_any(v0, r, pieces, n):
    return div_any(psi_series_any(v0 + r, pieces, n), psi_series_any(v0, pieces, n))


def biv_zero(ns, nu):
    return [[mp.mpf("0") for _ in range(nu)] for _ in range(ns)]


def biv_from_s(series, ns, nu):
    out = biv_zero(ns, nu)
    for i in range(min(ns, len(series))):
        out[i][0] = series[i]
    return out


def biv_from_u(series, ns, nu):
    out = biv_zero(ns, nu)
    for j in range(min(nu, len(series))):
        out[0][j] = series[j]
    return out


def biv_add(a, b):
    return [[scalar_add(a[i][j], b[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


def biv_scale(a, k):
    return [[scalar_mul(a[i][j], k) for j in range(len(a[0]))] for i in range(len(a))]


def biv_mul(a, b):
    ns = len(a)
    nu = len(a[0])
    out = biv_zero(ns, nu)
    for i in range(ns):
        for j in range(nu):
            out[i][j] = sum_any(
                scalar_mul(a[p][q], b[i - p][j - q])
                for p in range(i + 1)
                for q in range(j + 1)
            )
    return out


def biv_exp(a):
    ns = len(a)
    nu = len(a[0])
    out = biv_zero(ns, nu)
    out[0][0] = exp_any(a[0][0])
    # Recurrence from D exp(A)=exp(A) D A.  Works coefficientwise for the
    # bivariate ordinary Taylor coefficients used here.
    for total in range(1, ns + nu - 1):
        for i in range(ns):
            j = total - i
            if j < 0 or j >= nu:
                continue
            acc = None
            for p in range(i + 1):
                for q in range(j + 1):
                    if p == 0 and q == 0:
                        continue
                    if p == 0:
                        factor = q
                    elif q == 0:
                        factor = p
                    else:
                        factor = p
                    # Use the s-recursion when possible, otherwise u-recursion.
                    if i > 0 and p <= i and q <= j:
                        term = scalar_mul(scalar_mul(a[p][q], p), out[i - p][j - q])
                    elif i == 0 and q <= j:
                        term = scalar_mul(scalar_mul(a[p][q], q), out[i - p][j - q])
                    else:
                        continue
                    acc = term if acc is None else scalar_add(acc, term)
            out[i][j] = (mp.mpf("0") if acc is None else acc) / (i if i > 0 else j)
    return out


def weight_biv(center, omega, ns, nu):
    y = biv_zero(ns, nu)
    y[0][0] = center
    if ns > 1:
        y[1][0] = mp.mpf("0.5")
    if nu > 1:
        y[0][1] = mp.mpf("0.5")
    ep = biv_exp(biv_scale(y, omega))
    em = biv_exp(biv_scale(y, -omega))
    cosh = biv_scale(biv_add(ep, em), mp.mpf("0.5"))
    return biv_mul(y, cosh)


def h_derivatives_du_order_for_pieces(args, pieces, u, u_order, r_nodes, r_weights):
    ns = args.jet_order
    nu = u_order + 1
    s0 = mp.mpf(args.s0)
    omega = mp.mpf(args.omega)
    coeffs = [mp.mpf("0") for _ in range(ns)]
    for r, weight in zip(r_nodes, r_weights):
        ratio_s = biv_from_s(ratio_series_any(s0, r, pieces, ns), ns, nu)
        ratio_u = biv_from_u(ratio_series_any(u, r, pieces, nu), ns, nu)
        center = scalar_add(scalar_mul(u, mp.mpf("0.5")), r + s0 / 2)
        w = weight_biv(center, omega, ns, nu)
        term = biv_mul(biv_mul(w, ratio_s), ratio_u)
        for i in range(ns):
            coeffs[i] = scalar_add(coeffs[i], scalar_mul(term[i][u_order], weight))
    return [
        scalar_mul(coeffs[i], scalar_mul(mp.factorial(i), mp.factorial(u_order)))
        for i in range(ns)
    ]


def source_rows_du_order_for_pieces(args, polys, e_derivs, pieces, r_nodes, r_weights, u, u_order):
    h_derivs = h_derivatives_du_order_for_pieces(
        args,
        pieces,
        u,
        u_order,
        r_nodes,
        r_weights,
    )
    if hasattr(u, "_mpi_"):
        pstar = adjoint_source_value_interval(e_derivs, h_derivs, args.jet_order)
        brow = trace_green_concomitant_row_interval(e_derivs, h_derivs, args.jet_order)
        brow_poly = boundary_functional_row_interval(polys, mp.mpf(args.s0), brow)
        eval_poly = eval_functional_row_interval(polys, mp.mpf(args.s0), pstar)
        return [brow_poly, eval_poly]

    pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
    brow = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
    brow_poly = boundary_functional_row(polys, mp.mpf(args.s0), brow)
    eval_poly = eval_functional_row(polys, mp.mpf(args.s0), pstar)
    out = mp.matrix(2, len(polys))
    for j in range(len(polys)):
        out[0, j] = brow_poly[0, j]
        out[1, j] = eval_poly[0, j]
    return out


def source_row_scaled(args, base, e_derivs, pieces, r_nodes, r_weights, u, u_order):
    rows = source_rows_du_order_for_pieces(
        args,
        base["polys"],
        e_derivs,
        pieces,
        r_nodes,
        r_weights,
        u,
        u_order,
    )
    return rows * base["scaledModes"]


def interval_source_row_scaled(args, base, e_derivs, pieces, r_nodes, r_weights, u_iv, u_order):
    rows = source_rows_du_order_for_pieces(
        args,
        base["polys"],
        e_derivs,
        pieces,
        r_nodes,
        r_weights,
        u_iv,
        u_order,
    )
    out = []
    for row in rows:
        scaled = []
        for j in range(base["scaledModes"].cols):
            scaled.append(
                sum_any(row[k] * base["scaledModes"][k, j] for k in range(len(row)))
            )
        out.append(scaled)
    return out


def cross_interval(A, B):
    rows = len(A)
    cols = len(A[0])
    out = [[mp.mpf("0") for _ in range(cols)] for _ in range(cols)]
    for i in range(cols):
        for j in range(cols):
            out[i][j] = sum_any(A[k][i] * B[k][j] for k in range(rows))
    return out


def d2_interval_gram(E0, E1, E2):
    left = cross_interval(E2, E0)
    mid = cross_interval(E1, E1)
    right = cross_interval(E0, E2)
    out = [[mp.mpf("0") for _ in range(len(left))] for _ in range(len(left))]
    for i in range(len(out)):
        for j in range(len(out)):
            out[i][j] = left[i][j] + 2 * mid[i][j] + right[i][j]
    return out


def interval_envelope(args, base, e_derivs, pieces, r_nodes, r_weights):
    lo = mp.mpf(args.source_min)
    hi = mp.mpf(args.source_max)
    step = (hi - lo) / args.subintervals
    max_bound = mp.mpf("0")
    worst = None
    rows = []
    min_denom = mp.inf
    for k in range(args.subintervals):
        a = lo + k * step
        b = hi if k == args.subintervals - 1 else lo + (k + 1) * step
        u_iv = mp.iv.mpf([a, b])
        E0 = interval_source_row_scaled(args, base, e_derivs, pieces, r_nodes, r_weights, u_iv, 0)
        E1 = interval_source_row_scaled(args, base, e_derivs, pieces, r_nodes, r_weights, u_iv, 1)
        E2 = interval_source_row_scaled(args, base, e_derivs, pieces, r_nodes, r_weights, u_iv, 2)
        d2s = d2_interval_gram(E0, E1, E2)
        bound = interval_matrix_norm_bound(d2s)
        row = {
            "left": f(a),
            "right": f(b),
            "d2SourceGramBound": f(bound),
        }
        rows.append(row)
        if bound > max_bound:
            max_bound = bound
            worst = row
        # Cheap denominator lower diagnostic over endpoints and interval.
        for x in (a, b):
            val = psi_series_any(x, pieces, 1)[0]
            min_denom = min(min_denom, abs(val))
        if (k + 1) % max(1, args.subintervals // 8) == 0:
            print(f"  intervals {k+1:4d}/{args.subintervals} M2={fmt(max_bound, 8)}", flush=True)
    return max_bound, rows, worst, min_denom


def op_norm_sym(mat):
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def eigensystem_source_gram(source_gram):
    vals, vecs = mp.eigsy(sym(source_gram), eigvals_only=False)
    active_idx = [len(vals) - 2, len(vals) - 1]
    active = columns(vecs, active_idx)
    complement_idx = list(range(len(vals) - 2))
    complement = columns(vecs, complement_idx)
    gap = vals[-2] - vals[-3]
    return vals, active_idx, active, complement, gap


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
    parser.add_argument("--subintervals", type=int, default=32)
    parser.add_argument("--full-nmax", type=int, default=8)
    parser.add_argument("--source-order", type=int, default=10)
    parser.add_argument("--source-rmax", default="8")
    parser.add_argument("--tail-json", default="full_theta_interval_propagation.json")
    parser.add_argument(
        "--reuse-envelope-json",
        default=None,
        help="Reuse d_u^2(E_u^*E_u) interval envelope from a previous compatible run.",
    )
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
    parser.add_argument("--dps", type=int, default=55)
    parser.add_argument("--iv-dps", type=int, default=35)
    parser.add_argument("--matrix-order", type=int, default=40)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=16)
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
    parser.add_argument("--json-out", default="full_theta_source_quadrature_certificate.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps

    tail_data = json.loads(Path(args.tail_json).read_text(encoding="utf-8"))
    omitted_tail_bound = mp.mpf(tail_data["deltaSOperatorBound"])

    pieces8 = pieces_from_mode_weights(full_weights(args.full_nmax))
    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace/high block basis={args.basis} max_q={max_q}", flush=True)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)
    base_args = SimpleNamespace(**vars(args))
    base_args.source_grid = 9
    base = local_case(base_args, args.basis, e_derivs_at_s0)

    print(f"Building finite S8 source grid={args.source_grid}", flush=True)
    nodes, weights, e8, s8 = weighted_source_matrix_for_pieces(
        args,
        base,
        e_derivs_at_s0,
        pieces8,
    )
    s8_norm = op_norm_sym(s8)
    e8_norm = op_norm_rect(e8)
    svals, active_idx, active_basis, complement_basis, gap = eigensystem_source_gram(s8)
    lmat, finite_svals, rank_margin, response_norm = finite_riesz_margin(
        args,
        base,
        active_basis,
    )
    complement_response_norm = op_norm_rect(lmat * complement_basis)

    if args.reuse_envelope_json:
        env_data = json.loads(Path(args.reuse_envelope_json).read_text(encoding="utf-8"))
        m2 = mp.mpf(env_data["d2SourceGramIntervalEnvelope"])
        interval_rows = env_data.get("intervalRows", [])
        worst_interval = env_data.get("worstInterval", {})
        min_denom = mp.mpf(env_data.get("minRatioDenominatorDiagnostic", "nan"))
        print(
            f"Reusing interval d_u^2(E^*E) envelope from {args.reuse_envelope_json}",
            flush=True,
        )
    else:
        print("Bounding interval d_u^2(E^*E) envelope", flush=True)
        r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.source_rmax), args.source_order)
        m2, interval_rows, worst_interval, min_denom = interval_envelope(
            args,
            base,
            e_derivs_at_s0,
            pieces8,
            r_nodes,
            r_weights,
        )

    a = mp.mpf(args.source_min)
    b = mp.mpf(args.source_max)
    h = (b - a) / (args.source_grid - 1)
    trapezoid_factor = (b - a) * h**2 / 12
    quadrature_error = trapezoid_factor * m2
    total_error = quadrature_error + omitted_tail_bound
    alpha = 4 * total_error / gap if gap else mp.inf
    angle_factor = mp.sqrt(max(mp.mpf("0"), 1 - alpha**2)) if alpha < 1 else mp.mpf("0")
    lower_after_total = rank_margin * angle_factor - complement_response_norm * alpha
    pass_gap = bool(total_error < gap / 4)
    pass_rank = bool(pass_gap and lower_after_total > 0)

    print("Full-theta source quadrature certificate", flush=True)
    print(f"  M2={fmt(m2, 12)} factor={fmt(trapezoid_factor, 12)}", flush=True)
    print(f"  quadrature_error={fmt(quadrature_error, 12)}", flush=True)
    print(f"  omitted_tail_bound={fmt(omitted_tail_bound, 12)}", flush=True)
    print(f"  total_error={fmt(total_error, 12)} total/S8={fmt(total_error/s8_norm, 12)}", flush=True)
    print(f"  gap={fmt(gap, 12)} alpha={fmt(alpha, 12)}", flush=True)
    print(f"  ||L||={fmt(response_norm, 12)} ||L Pcomp||={fmt(complement_response_norm, 12)}", flush=True)
    print(f"  lower_after_total={fmt(lower_after_total, 12)} pass={pass_rank}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "sourceGrid": args.source_grid,
        "sourceNodesFirstLast": [f(nodes[0]), f(nodes[-1])],
        "sourceWeightsFirstLast": [f(weights[0]), f(weights[-1])],
        "subintervals": args.subintervals,
        "fullNmax": args.full_nmax,
        "sourceOrder": args.source_order,
        "sourceRmax": f(mp.mpf(args.source_rmax)),
        "tailJson": args.tail_json,
        "reuseEnvelopeJson": args.reuse_envelope_json,
        "sourceGramNormS8": f(s8_norm),
        "finiteE8OperatorNorm": f(e8_norm),
        "d2SourceGramIntervalEnvelope": f(m2),
        "trapezoidFactor": f(trapezoid_factor),
        "sourceQuadratureErrorBound": f(quadrature_error),
        "omittedTailDeltaSBound": f(omitted_tail_bound),
        "totalContinuumErrorBound": f(total_error),
        "totalContinuumErrorRelativeToS8": f(total_error / s8_norm if s8_norm else mp.inf),
        "activeEigenvaluesS8": [f(svals[i]) for i in active_idx],
        "complementTopEigenvalueS8": f(svals[-3]),
        "spectralGapToComplementS8": f(gap),
        "totalProjectorAlpha": f(alpha),
        "totalGapConditionPasses": pass_gap,
        "finiteActiveResponseSingularValuesS8": [f(x) for x in finite_svals],
        "finiteRankMarginS8": f(rank_margin),
        "responseOperatorNorm": f(response_norm),
        "complementResponseOperatorNorm": f(complement_response_norm),
        "angleFactorLowerBound": f(angle_factor),
        "lowerBoundAfterContinuumAndTail": f(lower_after_total),
        "fullPhiContinuumSourceNoncollapsePasses": pass_rank,
        "minRatioDenominatorDiagnostic": f(min_denom),
        "worstInterval": worst_interval,
        "intervalRows": interval_rows,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Continuum source-quadrature certificate for the literal full-theta "
            "source row.  The finite Psi_8 source Gram is promoted to the "
            "source window by an interval bound on d_u^2(E_u^*E_u); the omitted "
            "n>=9 interval-propagation bound is then added before applying the "
            "Riesz projector/rank margin."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
