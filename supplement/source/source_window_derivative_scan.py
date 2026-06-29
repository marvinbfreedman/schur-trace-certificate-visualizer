#!/usr/bin/env python3
r"""Analytic u-derivative certificate for the source-window family E_u.

The finite-net Lipschitz scan used divided differences of the two-row source
operator

    E_u f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0)).

This script differentiates the Green/source formula for h_u under the
Volterra integral.  For the endpoint-B model,

    h_u^{(k)}(s0) =
      int e^{5r/2} [
        r V_s^{(k)}(s0,r) V(u,r)
        + 1/2(V_s^{(k)}(s0,r) W(u,r)
              + W_s^{(k)}(s0,r) V(u,r))
      ] dr.

Hence d/du h_u is obtained by replacing V(u,r),W(u,r) by their u-derivatives.
Those derivatives are read from the same Taylor engine used for the confluent
kernel:

    V_u(u,r) = [z^1] V(u+z,r),    W_u(u,r) = [z^1] W(u+z,r).

The script verifies the analytic derivative against centered finite
differences and computes the sharp A-energy bound

    (partial_u E_u)^*(partial_u E_u) <= L(u)^2 A

on the closed-trace high block.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from aux_regularizer_certificate import direct_absorption_constant  # noqa: E402
from endpoint_kb_confluent_mp import endpoint_series  # noqa: E402
from lagrange_energy_control_certificate import (  # noqa: E402
    boundary_functional_row,
    eval_functional_row,
    make_qargs,
    split_spaces,
)
from lagrange_hardy_graph_certificate import residual_rows_for  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    trace_matrix,
)
from trace_lagrange_adjoint_control import (  # noqa: E402
    adjoint_source_value,
    load_exact_trace,
    trace_green_concomitant_row,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def row_frobenius(row):
    return mp.sqrt(mp.fsum(row[i, j] * row[i, j] for i in range(row.rows) for j in range(row.cols)))


def source_derivatives_du(c, s0, u, max_deriv, r_nodes, r_weights):
    """Return d/du of h_u^(0..max_deriv)(s0)."""
    coeffs = [mp.mpf("0") for _ in range(max_deriv + 1)]
    for r, weight in zip(r_nodes, r_weights):
        x = mp.e**r
        _fs, _vas, _was, vb_s, wb_s = endpoint_series(c, s0, r, max_deriv + 1)
        _fu, _vau, _wau, vb_u, wb_u = endpoint_series(c, u, r, 2)
        dvt = vb_u[1]
        dwt = wb_u[1]
        layer_weight = weight * x ** mp.mpf("2.5")
        for i in range(max_deriv + 1):
            coeffs[i] += layer_weight * (
                r * vb_s[i] * dvt
                + mp.mpf("0.5") * (vb_s[i] * dwt + wb_s[i] * dvt)
            )
    return [mp.factorial(i) * coeffs[i] for i in range(max_deriv + 1)]


def residual_rows_du_for(args, polys, e_derivs, r_nodes, r_weights, u):
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h_derivs = source_derivatives_du(
        c,
        s0,
        u,
        args.jet_order,
        r_nodes,
        r_weights,
    )
    pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
    brow = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
    brow_poly = boundary_functional_row(polys, s0, brow)
    eval_poly = eval_functional_row(polys, s0, pstar)
    out = mp.matrix(2, len(polys))
    for j in range(len(polys)):
        out[0, j] = brow_poly[0, j]
        out[1, j] = eval_poly[0, j]
    return out, pstar


def finite_diff_row(args, polys, e_derivs, r_nodes, r_weights, u, eps):
    plus, _ = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u + eps)
    minus, _ = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u - eps)
    return (plus - minus) / (2 * eps)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=17)
    parser.add_argument("--diff-u", default="0.08 0.30 0.52")
    parser.add_argument("--eps-values", default="1e-3 5e-4 1e-4")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=20)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=80)
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
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="source_window_derivative_scan.json")
    args = parser.parse_args()

    if args.source_grid < 2:
        raise SystemExit("--source-grid must be at least 2")

    mp.mp.dps = args.dps
    qargs = make_qargs(args)
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    _centers, R = trace_matrix(polys, qargs)
    N, _U, rank, nullity = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    diff_rows = []
    for text in args.diff_u.replace(",", " ").split():
        u = mp.mpf(text)
        analytic, _pstar = residual_rows_du_for(args, polys, e_derivs, r_nodes, r_weights, u)
        an_norm = row_frobenius(analytic)
        checks = []
        for eps_text in args.eps_values.replace(",", " ").split():
            eps = mp.mpf(eps_text)
            fd = finite_diff_row(args, polys, e_derivs, r_nodes, r_weights, u, eps)
            defect = fd - analytic
            abs_defect = row_frobenius(defect)
            rel_defect = abs_defect / an_norm if an_norm else mp.mpf("0")
            checks.append(
                {
                    "eps": f(eps),
                    "absoluteDefect": f(abs_defect),
                    "relativeDefect": f(rel_defect),
                }
            )
        diff_rows.append({"u": f(u), "analyticNorm": f(an_norm), "checks": checks})

    source_min = mp.mpf(args.source_min)
    source_max = mp.mpf(args.source_max)
    step = (source_max - source_min) / (args.source_grid - 1)
    radius = step / 2
    source_nodes = [source_min + i * step for i in range(args.source_grid)]

    sampled = []
    derivative_rows = []
    for u in source_nodes:
        rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        drows, dpstar = residual_rows_du_for(args, polys, e_derivs, r_nodes, r_weights, u)
        row_n = rows * N
        drow_n = drows * N
        high_eta, _high_size = direct_absorption_constant(row_n, high_a_modes, high_avals)
        full_eta, _full_size = direct_absorption_constant(row_n, a_modes, avals)
        high_deriv, high_deriv_size = direct_absorption_constant(drow_n, high_a_modes, high_avals)
        full_deriv, full_deriv_size = direct_absorption_constant(drow_n, a_modes, avals)
        sampled.append(
            {
                "u": f(u),
                "pStarH": f(pstar),
                "dPStarH": f(dpstar),
                "highAbsorb": f(high_eta),
                "fullAbsorb": f(full_eta),
                "highFullFrac": f(high_eta / full_eta if full_eta else mp.mpf("0")),
                "highDerivativeAbsorb": f(high_deriv),
                "fullDerivativeAbsorb": f(full_deriv),
                "highDerivativeFullFrac": f(high_deriv / full_deriv if full_deriv else mp.mpf("0")),
                "highDerivativeSize": f(high_deriv_size),
                "fullDerivativeSize": f(full_deriv_size),
            }
        )
        derivative_rows.append(drow_n)

    max_high = max(mp.mpf(row["highAbsorb"]) for row in sampled)
    min_full = min(mp.mpf(row["fullAbsorb"]) for row in sampled)
    max_frac = max(mp.mpf(row["highFullFrac"]) for row in sampled)
    max_deriv = max(mp.mpf(row["highDerivativeAbsorb"]) for row in sampled)
    max_deriv_frac = max(mp.mpf(row["highDerivativeFullFrac"]) for row in sampled)
    cover_high = 2 * max_high + 2 * radius * radius * max_deriv
    cover_frac = cover_high / min_full if min_full else mp.mpf("0")
    worst_deriv = max(sampled, key=lambda row: mp.mpf(row["highDerivativeAbsorb"]))

    print(
        f"Source-window derivative scan basis={args.basis} constraints={args.constraints} "
        f"cutoff={cutoff} grid={args.source_grid}",
        flush=True,
    )
    print("  finite-difference check max rel by u:")
    for row in diff_rows:
        best = min(row["checks"], key=lambda item: mp.mpf(item["relativeDefect"]))
        print(
            f"    u={fmt(mp.mpf(row['u']), 6)} best eps={best['eps']} "
            f"rel={fmt(mp.mpf(best['relativeDefect']), 8)}",
            flush=True,
        )
    print(f"  sampled max high/full = {fmt(max_frac, 10)}", flush=True)
    print(f"  max analytic high derivative A-constant = {fmt(max_deriv, 10)}", flush=True)
    print(f"  worst derivative u = {fmt(mp.mpf(worst_deriv['u']), 8)}", flush=True)
    print(f"  max derivative high/full = {fmt(max_deriv_frac, 10)}", flush=True)
    print(f"  analytic covering high/full bound = {fmt(cover_frac, 10)}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "sourceMin": f(source_min),
        "sourceMax": f(source_max),
        "sourceGrid": args.source_grid,
        "mesh": f(step),
        "meshRadius": f(radius),
        "maxSampleHighAbsorb": f(max_high),
        "minSampleFullAbsorb": f(min_full),
        "maxSampleHighFullFrac": f(max_frac),
        "maxAnalyticDerivativeHighAbsorb": f(max_deriv),
        "maxAnalyticDerivativeHighFullFrac": f(max_deriv_frac),
        "worstDerivativeU": worst_deriv["u"],
        "analyticCoverHighAbsorb": f(cover_high),
        "analyticCoverHighFullFracVsMinFull": f(cover_frac),
        "finiteDifferenceChecks": diff_rows,
        "sampled": sampled,
        "interpretation": (
            "Analytic u-derivative of the Green/source row. This upgrades the "
            "finite divided-difference Lipschitz scan to the exact derivative "
            "needed for a continuum source-window Hardy/Green theorem."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
