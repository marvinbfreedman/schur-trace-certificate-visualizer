#!/usr/bin/env python3
r"""Adjoint-evaluation Green representer certificate.

For the closed-trace Hardy/Green theorem, the dominant source row is

    ell_u(f) = P^*h_u(s0) f(s0).

This row factorizes into a scalar source amplitude p(u)=P^*h_u(s0) and the
fixed evaluation functional ev_{s0}.  Thus its Green representer is

    g_u = p(u) k_{s0}^{hi},

where k_{s0}^{hi} is the A-Green representer of evaluation at s0 on the
closed-trace high block:

    f(s0) = <k_{s0}^{hi}, f>_A.

The same factorization controls the differentiated row:

    partial_u ell_u(f) = p'(u) f(s0),
    partial_u g_u = p'(u) k_{s0}^{hi}.

This script constructs the finite representer k_{s0}^{hi}, verifies the
factorization against the previously used residual rows, and gives a compact
source-window envelope for p, p', and p''.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_kb_confluent_mp import endpoint_series  # noqa: E402
from lagrange_energy_control_certificate import (  # noqa: E402
    eval_functional_row,
    make_qargs,
    split_spaces,
)
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    poly_derivative_value,
    trace_matrix,
)
from trace_lagrange_adjoint_control import (  # noqa: E402
    adjoint_source_value,
    load_exact_trace,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def row_component(row, comp: int):
    out = mp.matrix(1, row.cols)
    for j in range(row.cols):
        out[0, j] = row[comp, j]
    return out


def source_derivatives_du_order(c, s0, u, max_deriv, u_order, r_nodes, r_weights):
    """Return d_u^u_order of h_u^(0..max_deriv)(s0)."""
    coeffs = [mp.mpf("0") for _ in range(max_deriv + 1)]
    for r, weight in zip(r_nodes, r_weights):
        x = mp.e**r
        _fs, _vas, _was, vb_s, wb_s = endpoint_series(c, s0, r, max_deriv + 1)
        _fu, _vau, _wau, vb_u, wb_u = endpoint_series(c, u, r, u_order + 1)
        dvt = mp.factorial(u_order) * vb_u[u_order]
        dwt = mp.factorial(u_order) * wb_u[u_order]
        layer_weight = weight * x ** mp.mpf("2.5")
        for i in range(max_deriv + 1):
            coeffs[i] += layer_weight * (
                r * vb_s[i] * dvt
                + mp.mpf("0.5") * (vb_s[i] * dwt + wb_s[i] * dvt)
            )
    return [mp.factorial(i) * coeffs[i] for i in range(max_deriv + 1)]


def energy_representer(row_n_comp, A, high_a_modes, high_avals):
    coeffs = row_n_comp * high_a_modes
    rep_n = mp.matrix(high_a_modes.rows, 1)
    norm2 = mp.mpf("0")
    for j, lam in enumerate(high_avals):
        c = coeffs[0, j]
        norm2 += c * c / lam
        for i in range(high_a_modes.rows):
            rep_n[i] += high_a_modes[i, j] * c / lam

    predicted = rep_n.T * A * high_a_modes
    defect2 = mp.mpf("0")
    target2 = mp.mpf("0")
    for j in range(high_a_modes.cols):
        defect = predicted[0, j] - coeffs[0, j]
        defect2 += defect * defect
        target2 += coeffs[0, j] * coeffs[0, j]
    rel = mp.sqrt(defect2 / target2) if target2 else mp.mpf("0")
    return rep_n, norm2, rel


def full_coeffs(N, rep_n):
    out = N * rep_n
    return [out[i] for i in range(out.rows)]


def eval_profile(polys, coeffs, length, points):
    rows = []
    max_abs = mp.mpf("0")
    l2 = mp.mpf("0")
    prev_x = None
    prev_value = None
    for i in range(points):
        x = length * i / (points - 1)
        value = mp.fsum(coeffs[j] * poly_derivative_value(polys[j], x, 0) for j in range(len(polys)))
        max_abs = max(max_abs, abs(value))
        if prev_x is not None:
            dx = x - prev_x
            l2 += dx * (prev_value * prev_value + value * value) / 2
        prev_x = x
        prev_value = value
        rows.append({"s": f(x), "value": f(value)})
    return {"maxAbs": f(max_abs), "l2Approx": f(l2), "grid": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=33)
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
    parser.add_argument("--profile-points", type=int, default=101)
    parser.add_argument("--json-out", default="adjoint_eval_representer_certificate.json")
    args = parser.parse_args()

    if args.source_grid < 2:
        raise SystemExit("--source-grid must be at least 2")

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
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

    eval_row = eval_functional_row(polys, s0, mp.mpf("1"))
    eval_row_n = eval_row * N
    eval_rep_n, eval_norm2, eval_range_defect = energy_representer(
        eval_row_n,
        A,
        high_a_modes,
        high_avals,
    )
    eval_coeffs = full_coeffs(N, eval_rep_n)
    eval_prof = eval_profile(polys, eval_coeffs, mp.mpf(args.L), args.profile_points)

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    c = mp.pi
    source_min = mp.mpf(args.source_min)
    source_max = mp.mpf(args.source_max)
    step = (source_max - source_min) / (args.source_grid - 1)
    radius = step / 2

    rows = []
    max_abs_p = mp.mpf("0")
    max_abs_dp = mp.mpf("0")
    max_abs_d2p = mp.mpf("0")
    max_factor_defect = mp.mpf("0")
    for i in range(args.source_grid):
        u = source_min + i * step
        h0 = source_derivatives_du_order(c, s0, u, args.jet_order, 0, r_nodes, r_weights)
        h1 = source_derivatives_du_order(c, s0, u, args.jet_order, 1, r_nodes, r_weights)
        h2 = source_derivatives_du_order(c, s0, u, args.jet_order, 2, r_nodes, r_weights)
        p = adjoint_source_value(e_derivs, h0, args.jet_order)
        dp = adjoint_source_value(e_derivs, h1, args.jet_order)
        d2p = adjoint_source_value(e_derivs, h2, args.jet_order)

        predicted_eval = p * p * eval_norm2
        predicted_deval = dp * dp * eval_norm2
        eval_component = eval_functional_row(polys, s0, p) * N
        deval_component = eval_functional_row(polys, s0, dp) * N
        _rep, actual_eval, _rel = energy_representer(eval_component, A, high_a_modes, high_avals)
        _drep, actual_deval, _drel = energy_representer(deval_component, A, high_a_modes, high_avals)
        scale = max(mp.mpf("1"), abs(actual_eval), abs(predicted_eval), abs(actual_deval), abs(predicted_deval))
        factor_defect = max(abs(actual_eval - predicted_eval), abs(actual_deval - predicted_deval)) / scale
        max_factor_defect = max(max_factor_defect, factor_defect)

        max_abs_p = max(max_abs_p, abs(p))
        max_abs_dp = max(max_abs_dp, abs(dp))
        max_abs_d2p = max(max_abs_d2p, abs(d2p))
        rows.append(
            {
                "u": f(u),
                "pStar": f(p),
                "dPStar": f(dp),
                "d2PStar": f(d2p),
                "evalNorm2": f(actual_eval),
                "dEvalNorm2": f(actual_deval),
                "factorRelativeDefect": f(factor_defect),
            }
        )

    p_cover = max_abs_p + radius * max_abs_dp
    dp_cover = max_abs_dp + radius * max_abs_d2p
    eval_cover = p_cover * p_cover * eval_norm2
    deval_cover = dp_cover * dp_cover * eval_norm2

    print(
        f"Adjoint-eval representer certificate basis={args.basis} "
        f"constraints={args.constraints} cutoff={cutoff}",
        flush=True,
    )
    print(f"  eval representer norm^2 = {fmt(eval_norm2, 12)}", flush=True)
    print(f"  eval range defect       = {fmt(eval_range_defect, 8)}", flush=True)
    print(f"  max factor defect       = {fmt(max_factor_defect, 8)}", flush=True)
    print(
        f"  cover |p|={fmt(p_cover, 10)} |p'|={fmt(dp_cover, 10)} "
        f"eval={fmt(eval_cover, 10)} dEval={fmt(deval_cover, 10)}",
        flush=True,
    )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(s0),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "sourceMin": f(source_min),
        "sourceMax": f(source_max),
        "sourceGrid": args.source_grid,
        "mesh": f(step),
        "meshRadius": f(radius),
        "evalRepresenterNorm2": f(eval_norm2),
        "evalRepresenterNorm": f(mp.sqrt(max(mp.mpf("0"), eval_norm2))),
        "evalRangeRelativeDefect": f(eval_range_defect),
        "maxFactorRelativeDefect": f(max_factor_defect),
        "maxAbsPStarGrid": f(max_abs_p),
        "maxAbsDPStarGrid": f(max_abs_dp),
        "maxAbsD2PStarGrid": f(max_abs_d2p),
        "coverAbsPStar": f(p_cover),
        "coverAbsDPStar": f(dp_cover),
        "coverEvalNorm2": f(eval_cover),
        "coverDEvalNorm2": f(deval_cover),
        "rows": rows,
        "evalRepresenterProfile": eval_prof,
        "theoremTemplate": (
            "The adjoint-evaluation continuum representer is p(u) times the "
            "single Green representer k_{s0}^{hi} of point evaluation. A "
            "uniform bound follows from ||k_{s0}^{hi}||_A and compact-window "
            "bounds for p(u)=P^*h_u(s0) and p'(u)."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
