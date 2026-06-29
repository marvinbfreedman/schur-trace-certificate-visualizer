#!/usr/bin/env python3
r"""Complex-ball Bernstein-ellipse certificate for the source quadrature tail.

This is the rigorous version of ``source_quadrature_bernstein_bound.py``.  It
covers the Bernstein ellipse

    z(theta) = (rho + rho^{-1}) cos(theta)/2
             + i (rho - rho^{-1}) sin(theta)/2

by interval arcs, maps them to ``u=.30+.22 z``, and evaluates the explicit
Green-source formula with complex interval arithmetic.  The resulting bound

    M_rho >= sup_{z in E_rho} ||d_u^2(E(u)^*E(u))||

replaces the sampled ellipse proxy in the Bernstein tail estimate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_quadrature_bernstein_bound import (  # noqa: E402
    formal_pole_rho,
    load_finite_envelope,
)
from source_quadrature_error_bound import (  # noqa: E402
    op_norm_symmetric,
    trapezoid_source_gram,
)
from source_quadrature_interval_envelope import (  # noqa: E402
    d2_source_gram_interval,
    source_row_interval,
)


def iv_upper(x):
    return mp.mpf(x.b)


def iv_lower(x):
    return mp.mpf(x.a)


def complex_interval_abs_sup(x):
    if hasattr(x, "_mpci_") or type(x).__name__.endswith("ivmpc"):
        return iv_upper(abs(x))
    if hasattr(x, "_mpi_") or type(x).__name__.endswith("ivmpf"):
        return max(abs(iv_lower(x)), abs(iv_upper(x)))
    return abs(x)


def complex_interval_matrix_norm_bound(mat):
    row_bound = mp.mpf("0")
    frob2 = mp.mpf("0")
    for row in mat:
        rowsum = mp.mpf("0")
        for value in row:
            bound = complex_interval_abs_sup(value)
            rowsum += bound
            frob2 += bound * bound
        row_bound = max(row_bound, rowsum)
    frob = mp.sqrt(frob2)
    return min(row_bound, frob), row_bound, frob


def ellipse_arc_box(theta_left, theta_right, rho):
    theta = mp.iv.mpf([theta_left, theta_right])
    major = mp.iv.mpf([(rho + 1 / rho) / 2, (rho + 1 / rho) / 2])
    minor = mp.iv.mpf([(rho - 1 / rho) / 2, (rho - 1 / rho) / 2])
    real = major * mp.iv.cos(theta)
    imag = minor * mp.iv.sin(theta)
    return mp.iv.mpc(real, imag)


def interval_center_width(x):
    lo = iv_lower(x)
    hi = iv_upper(x)
    return (lo + hi) / 2, hi - lo


def arc_bound(args, base, e_derivs, r_nodes, r_weights, center, radius, rho, left, right):
    z_box = ellipse_arc_box(left, right, rho)
    center_iv = mp.iv.mpf([center, center])
    radius_iv = mp.iv.mpf([radius, radius])
    u_box = center_iv + radius_iv * z_box
    E0 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_box, 0)
    E1 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_box, 1)
    E2 = source_row_interval(args, base, e_derivs, r_nodes, r_weights, u_box, 2)
    d2s = d2_source_gram_interval(E0, E1, E2)
    bound, row_bound, frob_bound = complex_interval_matrix_norm_bound(d2s)
    _uc, real_width = interval_center_width(u_box.real)
    _vc, imag_width = interval_center_width(u_box.imag)
    return {
        "thetaLeft": f(left),
        "thetaRight": f(right),
        "uRealMin": f(iv_lower(u_box.real)),
        "uRealMax": f(iv_upper(u_box.real)),
        "uImagMin": f(iv_lower(u_box.imag)),
        "uImagMax": f(iv_upper(u_box.imag)),
        "uRealWidth": f(real_width),
        "uImagWidth": f(imag_width),
        "d2SourceGramBound": f(bound),
        "rowSumBound": f(row_bound),
        "frobeniusBound": f(frob_bound),
    }, bound, row_bound, frob_bound


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--finite-envelope-json", default="source_quadrature_chebyshev_ball_d32_g129.json")
    parser.add_argument("--rho", default="2.0")
    parser.add_argument("--arcs", type=int, default=64)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=129)
    parser.add_argument("--cheb-degree", type=int, default=32)
    parser.add_argument("--source-gram-norm", default="0")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
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
    parser.add_argument("--iv-dps", type=int, default=50)
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
    parser.add_argument("--json-out", default="source_quadrature_complex_ball_bernstein_rho2.json")
    args = parser.parse_args()

    if args.arcs < 4:
        raise SystemExit("--arcs must be at least 4")

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps

    finite_envelope, finite_key, finite_data = load_finite_envelope(args.finite_envelope_json)
    source_norm = mp.mpf(args.source_gram_norm)
    if source_norm <= 0 and "sourceGramNorm" in finite_data:
        source_norm = mp.mpf(finite_data["sourceGramNorm"])
    if source_norm <= 0:
        qargs = SimpleNamespace(**vars(args))
        qargs.source_grid = args.source_grid
        max_q_for_norm = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
        _vals_norm, e_derivs_norm, _lam_derivs_norm = exact_trace_derivatives(args, max_q_for_norm)
        base_norm = local_case(args, args.basis, e_derivs_norm)
        _nodes, _weights, s_quad = trapezoid_source_gram(qargs, base_norm, e_derivs_norm, args.source_grid)
        source_norm = op_norm_symmetric(s_quad)

    a = mp.mpf(args.source_min)
    b = mp.mpf(args.source_max)
    center = (a + b) / 2
    radius = (b - a) / 2
    h = (b - a) / (args.source_grid - 1)
    trapezoid_factor = (b - a) * h**2 / 12
    rho = mp.mpf(args.rho)
    tail_factor = 2 * rho ** (-(args.cheb_degree + 1)) / (1 - 1 / rho)
    target = mp.mpf(args.target_s_rel)
    allowed_total_envelope = target * source_norm / trapezoid_factor
    allowed_mrho = (allowed_total_envelope - finite_envelope) / tail_factor
    x0, rho_formula_pole = formal_pole_rho(a, b, mp.mpf("0"))

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace data max_q={max_q} basis={args.basis}", flush=True)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print("Preparing local high block", flush=True)
    base = local_case(args, args.basis, e_derivs)
    print("Preparing endpoint quadrature", flush=True)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    print(
        f"Complex-ball Bernstein source envelope basis={args.basis} "
        f"rho={rho} arcs={args.arcs}",
        flush=True,
    )
    rows = []
    max_bound = mp.mpf("0")
    max_row_bound = mp.mpf("0")
    max_frob_bound = mp.mpf("0")
    worst = None
    step = 2 * mp.pi / args.arcs
    for k in range(args.arcs):
        left = k * step
        right = (k + 1) * step
        row, bound, row_bound, frob_bound = arc_bound(
            args,
            base,
            e_derivs,
            r_nodes,
            r_weights,
            center,
            radius,
            rho,
            left,
            right,
        )
        rows.append(row)
        if bound > max_bound:
            max_bound = bound
            worst = row
        max_row_bound = max(max_row_bound, row_bound)
        max_frob_bound = max(max_frob_bound, frob_bound)
        if (k + 1) % max(1, args.arcs // 8) == 0:
            print(f"  arcs {k+1:4d}/{args.arcs} M={fmt(max_bound, 8)}", flush=True)

    bernstein_tail = max_bound * tail_factor
    total_envelope = finite_envelope + bernstein_tail
    error_bound = trapezoid_factor * total_envelope
    rel = error_bound / source_norm
    pass_target = bool(rel < target)
    pass_allowed = bool(max_bound < allowed_mrho)

    print(f"  finite_envelope={fmt(finite_envelope, 12)} ({finite_key})", flush=True)
    print(f"  complex_ball_Mrho={fmt(max_bound, 12)} allowed={fmt(allowed_mrho, 12)}", flush=True)
    print(f"  bernstein_tail={fmt(bernstein_tail, 12)}", flush=True)
    print(f"  rel_error={fmt(rel, 12)} target={target} pass={pass_target}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "sourceMin": f(a),
        "sourceMax": f(b),
        "sourceCenter": f(center),
        "sourceRadius": f(radius),
        "sourceGrid": args.source_grid,
        "chebDegree": args.cheb_degree,
        "rho": f(rho),
        "arcs": args.arcs,
        "mpDps": args.dps,
        "ivDps": args.iv_dps,
        "mappedFormulaPoleAtU0": f(x0),
        "rhoToFormulaPoleU0": f(rho_formula_pole),
        "trapezoidFactor": f(trapezoid_factor),
        "sourceGramNorm": f(source_norm),
        "finiteEnvelopeSource": args.finite_envelope_json,
        "finiteEnvelopeKey": finite_key,
        "finiteCoefficientEnvelope": f(finite_envelope),
        "complexBallMrho": f(max_bound),
        "complexBallRowMrho": f(max_row_bound),
        "complexBallFrobeniusMrho": f(max_frob_bound),
        "allowedMrhoForTarget": f(allowed_mrho),
        "complexBallMrhoPassesAllowed": pass_allowed,
        "bernsteinTailFactor": f(tail_factor),
        "bernsteinTailEnvelope": f(bernstein_tail),
        "bernsteinTotalEnvelope": f(total_envelope),
        "sourceGramErrorBound": f(error_bound),
        "sourceGramRelativeErrorBound": f(rel),
        "sourceGramTargetRelative": f(target),
        "sourceGramPassesTarget": pass_target,
        "worstArc": worst,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "arcRows": rows,
        "interpretation": (
            "Complex-interval cover of the Bernstein ellipse.  This replaces "
            "the sampled M_rho proxy by a true interval supremum over the "
            "explicit finite Green-source formula."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
