#!/usr/bin/env python3
r"""Interval certificate for the source-quadrature derivative envelope.

This replaces the sampled envelope in ``source_quadrature_error_bound.py`` by
an interval-arithmetic bound in the source variable ``u``.  The r-integral and
Galerkin matrices are still the finite certified objects used throughout the
current notes; the new rigor is the supremum over ``u in [0.08,0.52]``.

For the A-normalized high block, let

    E(u): H -> R^2,
    S = int E(u)^* E(u) du.

The composite trapezoid theorem gives

    ||S-S_h|| <= (b-a) h^2 / 12 sup_u ||d_u^2(E(u)^*E(u))||.

The script encloses the derivative matrix on a partition of the source window
and uses an absolute row-sum/Frobenius bound for the 2x2 high-block matrix.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders, submatrix_rows, sym  # noqa: E402
from endpoint_kb_confluent_mp import zero  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_quadrature_error_bound import op_norm_symmetric, trapezoid_source_gram  # noqa: E402
from trace_active_derivative_rank import derivative_matrix  # noqa: E402


def exp_any(x):
    return mp.iv.exp(x) if hasattr(x, "_mpi_") else mp.e**x


def sum_any(values):
    total = None
    for value in values:
        total = value if total is None else total + value
    return mp.mpf("0") if total is None else total


def scalar_sub_any(x, y):
    if hasattr(x, "_mpi_"):
        return x - y
    if hasattr(y, "_mpi_"):
        return (-y) + x
    return x - y


def scalar_mul_any(x, y):
    if hasattr(x, "_mpi_"):
        return x * y
    if hasattr(y, "_mpi_"):
        return y * x
    return x * y


def sub_poly_any(a, b):
    return [scalar_sub_any(x, y) for x, y in zip(a, b)]


def mul_any(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = sum_any(scalar_mul_any(a[j], b[i - j]) for j in range(i + 1))
    return out


def div_any(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = scalar_sub_any(
            a[i],
            sum_any(b[j] * out[i - j] for j in range(1, i + 1)),
        ) / b[0]
    return out


def exp_series_any(a):
    n = len(a)
    out = zero(n)
    out[0] = exp_any(a[0])
    for i in range(1, n):
        out[i] = sum_any(a[k] * (k * out[i - k]) for k in range(1, i + 1)) / i
    return out


def lambda_series_any(c, s0, n):
    lam0 = exp_any(s0) * c
    return [lam0 / mp.factorial(i) for i in range(n)]


def scale_any(a, k):
    return [x * k for x in a]


def endpoint_series_any(c, s0, r, n):
    x = mp.e**r
    tau = x - 1
    lam = lambda_series_any(c, s0, n)
    denom = lam[:]
    denom[0] -= c

    exp_lam = exp_series_any(scale_any(lam, -tau))
    exp_c = mp.e ** (-c * tau)

    one = zero(n)
    one[0] = 1
    sser = zero(n)
    sser[0] = s0
    if n > 1:
        sser[1] = 1

    f_num = exp_lam[:]
    f_num[0] -= exp_c
    f = div_any(f_num, denom)

    gb_pref = sub_poly_any(one, scale_any(lam, x))
    gb = mul_any(gb_pref, exp_lam)
    gb_c = (1 - c * x) * exp_c
    vb_num = gb[:]
    vb_num[0] -= gb_c
    vb = div_any(vb_num, denom)
    wb = div_any(mul_any(sser, gb), denom)

    ga_pref = sub_poly_any(scale_any(one, mp.mpf("1.5")), scale_any(lam, x))
    ga = mul_any(ga_pref, exp_lam)
    ga_c = (mp.mpf("1.5") - c * x) * exp_c
    va_num = ga[:]
    va_num[0] -= ga_c
    va = div_any(va_num, denom)
    wa = div_any(mul_any(sser, ga), denom)

    return f, va, wa, vb, wb


def source_derivatives_du_order_interval(c, s0, u, max_deriv, u_order, r_nodes, r_weights):
    coeffs = [None for _ in range(max_deriv + 1)]
    for r, weight in zip(r_nodes, r_weights):
        x = mp.e**r
        _fs, _vas, _was, vb_s, wb_s = endpoint_series_any(c, s0, r, max_deriv + 1)
        _fu, _vau, _wau, vb_u, wb_u = endpoint_series_any(c, u, r, u_order + 1)
        dvt = vb_u[u_order] * mp.factorial(u_order)
        dwt = wb_u[u_order] * mp.factorial(u_order)
        layer_weight = weight * x ** mp.mpf("2.5")
        for i in range(max_deriv + 1):
            expr = (
                dvt * (r * vb_s[i])
                + (dwt * vb_s[i] + dvt * wb_s[i]) * mp.mpf("0.5")
            )
            term = expr * layer_weight
            coeffs[i] = term if coeffs[i] is None else coeffs[i] + term
    return [
        (mp.mpf("0") if coeffs[i] is None else coeffs[i]) * mp.factorial(i)
        for i in range(max_deriv + 1)
    ]


def product_derivative_interval(e_derivs, h_derivs, k: int, n: int):
    return sum_any(
        h_derivs[n - m]
        * (mp.binomial(n, m) * e_derivs[m][k] / mp.factorial(k))
        for m in range(n + 1)
    )


def adjoint_source_value_interval(e_derivs, h_derivs, jet_order: int):
    return sum_any(
        product_derivative_interval(e_derivs, h_derivs, k, k) * ((-1) ** k)
        for k in range(jet_order)
    )


def trace_green_concomitant_row_interval(e_derivs, h_derivs, jet_order: int):
    order = jet_order - 1
    row = [None for _ in range(order)]
    for k in range(1, jet_order):
        for j in range(k):
            n = k - 1 - j
            term = product_derivative_interval(e_derivs, h_derivs, k, n) * ((-1) ** n)
            row[j] = term if row[j] is None else row[j] + term
    return [mp.mpf("0") if value is None else value for value in row]


def boundary_functional_row_interval(polys, s0, brow):
    row = []
    for i, poly in enumerate(polys):
        total = None
        for j in range(len(brow)):
            term = brow[j] * poly_derivative_value_local(poly, s0, j)
            total = term if total is None else total + term
        row.append(mp.mpf("0") if total is None else total)
    return row


def poly_derivative_value_local(poly, x, deriv):
    total = mp.mpf("0")
    for power in range(deriv, len(poly)):
        coeff = poly[power] * mp.factorial(power) / mp.factorial(power - deriv)
        total += coeff * x ** (power - deriv)
    return total


def eval_functional_row_interval(polys, s0, scale_value):
    row = []
    for i, poly in enumerate(polys):
        row.append(scale_value * poly_derivative_value_local(poly, s0, 0))
    return row


def residual_rows_du_order_interval(args, polys, e_derivs, r_nodes, r_weights, u, u_order):
    h_derivs = source_derivatives_du_order_interval(
        mp.pi,
        mp.mpf(args.s0),
        u,
        args.jet_order,
        u_order,
        r_nodes,
        r_weights,
    )
    pstar = adjoint_source_value_interval(e_derivs, h_derivs, args.jet_order)
    brow = trace_green_concomitant_row_interval(e_derivs, h_derivs, args.jet_order)
    brow_poly = boundary_functional_row_interval(polys, mp.mpf(args.s0), brow)
    eval_poly = eval_functional_row_interval(polys, mp.mpf(args.s0), pstar)
    return [brow_poly, eval_poly]


def interval_abs_sup(x):
    if hasattr(x, "_mpi_"):
        return max(abs(mp.mpf(x.a)), abs(mp.mpf(x.b)))
    return abs(x)


def interval_matrix_norm_bound(mat):
    row_bound = mp.mpf("0")
    frob2 = mp.mpf("0")
    for i in range(len(mat)):
        rowsum = mp.mpf("0")
        for j in range(len(mat[i])):
            a = interval_abs_sup(mat[i][j])
            rowsum += a
            frob2 += a * a
        row_bound = max(row_bound, rowsum)
    return min(row_bound, mp.sqrt(frob2))


def source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_iv, u_order):
    rows = residual_rows_du_order_interval(
        args,
        base["polys"],
        e_derivs,
        r_nodes,
        r_weights,
        u_iv,
        u_order,
    )
    out = []
    for i in range(len(rows)):
        out_row = []
        for j in range(base["scaledModes"].cols):
            out_row.append(
                sum_any(rows[i][k] * base["scaledModes"][k, j] for k in range(len(rows[i])))
            )
        out.append(out_row)
    return out


def cross_gram_interval(A, B):
    a_rows = len(A)
    a_cols = len(A[0])
    b_cols = len(B[0])
    out = [[mp.mpf("0") for _ in range(b_cols)] for _ in range(a_cols)]
    for i in range(a_cols):
        for j in range(b_cols):
            out[i][j] = sum_any(A[k][i] * B[k][j] for k in range(a_rows))
    return out


def d2_source_gram_interval(E0, E1, E2):
    left = cross_gram_interval(E2, E0)
    mid = cross_gram_interval(E1, E1)
    right = cross_gram_interval(E0, E2)
    out = [[mp.mpf("0") for _ in range(len(left[0]))] for _ in range(len(left))]
    for i in range(len(out)):
        for j in range(len(out[i])):
            out[i][j] = left[i][j] + mid[i][j] * 2 + right[i][j]
    return out


def interval_row(args, base, e_derivs, r_nodes, r_weights, lo, hi):
    u_iv = mp.iv.mpf([lo, hi])
    E0 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_iv, 0)
    E1 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_iv, 1)
    E2 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_iv, 2)
    d2s = d2_source_gram_interval(E0, E1, E2)
    bound = interval_matrix_norm_bound(d2s)
    return {
        "left": f(lo),
        "right": f(hi),
        "width": f(hi - lo),
        "d2SourceGramBound": f(bound),
    }, bound


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=65)
    parser.add_argument("--subintervals", type=int, default=64)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
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
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--iv-dps", type=int, default=30)
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
    parser.add_argument("--target-s-rel", default="5.03e-6")
    parser.add_argument("--json-out", default="source_quadrature_interval_envelope.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs)

    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    qargs = SimpleNamespace(**vars(args))
    qargs.source_grid = args.source_grid
    _nodes, _weights, s_quad = trapezoid_source_gram(qargs, base, e_derivs, args.source_grid)
    s_norm = op_norm_symmetric(s_quad)

    lo = mp.mpf(args.source_min)
    hi = mp.mpf(args.source_max)
    h = (hi - lo) / (args.source_grid - 1)
    factor = (hi - lo) * h**2 / 12
    step = (hi - lo) / args.subintervals

    rows = []
    max_bound = mp.mpf("0")
    worst = None
    print(
        f"Interval source envelope basis={args.basis} source_grid={args.source_grid} "
        f"subintervals={args.subintervals}",
        flush=True,
    )
    for k in range(args.subintervals):
        a = lo + k * step
        b = hi if k == args.subintervals - 1 else lo + (k + 1) * step
        row, bound = interval_row(args, base, e_derivs, r_nodes, r_weights, a, b)
        rows.append(row)
        if bound > max_bound:
            max_bound = bound
            worst = row
        if (k + 1) % max(1, args.subintervals // 8) == 0:
            print(f"  {k+1:4d}/{args.subintervals} max={fmt(max_bound, 8)}", flush=True)

    error_bound = factor * max_bound
    rel = error_bound / s_norm if s_norm else mp.inf
    target = mp.mpf(args.target_s_rel)
    print(f"  envelope={fmt(max_bound, 12)}", flush=True)
    print(f"  ||S_h||={fmt(s_norm, 12)} rel_error={fmt(rel, 12)} target={target}", flush=True)
    print(f"  pass={rel < target}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "sourceMin": f(lo),
        "sourceMax": f(hi),
        "sourceGrid": args.source_grid,
        "subintervals": args.subintervals,
        "ivDps": args.iv_dps,
        "trapezoidFactor": f(factor),
        "sourceGramNorm": f(s_norm),
        "intervalSecondDerivativeEnvelope": f(max_bound),
        "sourceGramErrorBound": f(error_bound),
        "sourceGramRelativeErrorBound": f(rel),
        "sourceGramTargetRelative": f(target),
        "sourceGramPassesTarget": bool(rel < target),
        "worstInterval": worst,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "intervalRows": rows,
        "interpretation": (
            "Interval-arithmetic supremum bound for d_u^2(E^*E) on the source "
            "window.  This replaces the sampled derivative envelope in the "
            "trapezoid error theorem for S=E^*E."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
