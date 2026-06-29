#!/usr/bin/env python3
r"""Deterministic trapezoid-rule bounds for the weighted source kernels.

On the fixed basis-22 high block let

    E(u): H -> R^2

be the two-row source operator, with ``H`` A-normalized.  The continuum source
Gram operator is

    S = int_a^b E(u)^T E(u) du.

This is the right finite-dimensional representative of ``T=EE^*``: the
non-zero spectra agree.  Likewise ``B B^* = L S L^*`` for ``B=L E^*``.

For the composite trapezoid rule with mesh ``h``,

    ||S-S_h|| <= (b-a) h^2 / 12 sup_u ||d_u^2(E^T E)||.

The same bound applied after sandwiching by ``L`` gives a deterministic
quadrature estimate for ``B B^*``.  This script computes the analytic
derivatives of ``E(u)`` and evaluates the corresponding derivative envelopes.

The derivative envelope is sampled.  The displayed theorem is exact once the
reported envelope bound has been justified analytically, for example by the
same endpoint Green formula plus interval arithmetic or a Chebyshev tail bound.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_eval_representer_certificate import source_derivatives_du_order  # noqa: E402
from determinant_gap_bound_diagnostic import parse_orders, submatrix_rows, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints, stack_source_rows  # noqa: E402
from lagrange_energy_control_certificate import boundary_functional_row, eval_functional_row  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_side_quadrature_refinement import weights_for_nodes  # noqa: E402
from trace_active_derivative_rank import derivative_matrix  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value, trace_green_concomitant_row  # noqa: E402


def op_norm_symmetric(mat):
    if mat.rows == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def op_norm_rect(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat * mat.T), eigvals_only=True)
    return mp.sqrt(max([mp.mpf("0"), *vals]))


def residual_rows_du_order_for(args, polys, e_derivs, r_nodes, r_weights, u, u_order):
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h_derivs = source_derivatives_du_order(
        c,
        s0,
        u,
        args.jet_order,
        u_order,
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
    return out


def source_nodes_for(args, grid):
    lo = mp.mpf(args.source_min)
    hi = mp.mpf(args.source_max)
    if grid == 1:
        return [(lo + hi) / 2]
    step = (hi - lo) / (grid - 1)
    return [lo + i * step for i in range(grid)]


def source_row(args, base, e_derivs, r_nodes, r_weights, u, u_order):
    rows = residual_rows_du_order_for(args, base["polys"], e_derivs, r_nodes, r_weights, u, u_order)
    return rows * base["scaledModes"]


def second_derivative_integrand(E0, E1, E2):
    return E2.T * E0 + 2 * (E1.T * E1) + E0.T * E2


def fourth_derivative_integrand(E0, E1, E2, E3, E4):
    return (
        E4.T * E0
        + 4 * (E3.T * E1)
        + 6 * (E2.T * E2)
        + 4 * (E1.T * E3)
        + E0.T * E4
    )


def trapezoid_source_gram(args, base, e_derivs, grid):
    nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs)
    if len(nodes) != grid:
        # stack_source_rows reads args.source_grid; keep this explicit to avoid
        # accidentally mixing grids.
        raise ValueError(f"expected {grid} nodes, got {len(nodes)}")
    weights = weights_for_nodes(nodes)
    source_on_scaled = source_rows * base["scaledModes"]
    gram = mp.matrix(source_on_scaled.cols)
    for j, weight in enumerate(weights):
        block = mp.matrix(2, source_on_scaled.cols)
        for comp in range(2):
            for col in range(source_on_scaled.cols):
                block[comp, col] = source_on_scaled[2 * j + comp, col]
        gram += weight * block.T * block
    return nodes, weights, sym(gram)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=17)
    parser.add_argument("--envelope-grid", type=int, default=65)
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
    parser.add_argument("--target-b-rel", default="2.68e-4")
    parser.add_argument("--json-out", default="source_quadrature_error_bound.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs)

    point_args = SimpleNamespace(**vars(args))
    point_args.s0 = args.point
    _vals, _e_derivs_point, _lam_derivs = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(point_args, base, base["scaledModes"], mp.mpf(args.point), max(args.orders))
    lmat = submatrix_rows(full_response, args.orders)

    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    qargs = SimpleNamespace(**vars(args))
    qargs.source_grid = args.source_grid
    nodes, weights, s_quad = trapezoid_source_gram(qargs, base, e_derivs, args.source_grid)
    bgram_quad = sym(lmat * s_quad * lmat.T)

    a = mp.mpf(args.source_min)
    b = mp.mpf(args.source_max)
    h = (b - a) / (args.source_grid - 1)
    envelope_nodes = source_nodes_for(args, args.envelope_grid)

    max_s2 = mp.mpf("0")
    max_b2 = mp.mpf("0")
    max_s4 = mp.mpf("0")
    max_b4 = mp.mpf("0")
    max_e = mp.mpf("0")
    max_e1 = mp.mpf("0")
    max_e2 = mp.mpf("0")
    rows = []
    for u in envelope_nodes:
        E0 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 0)
        E1 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 1)
        E2 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 2)
        E3 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 3)
        E4 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 4)
        d2s = sym(second_derivative_integrand(E0, E1, E2))
        d2b = sym(lmat * d2s * lmat.T)
        d4s = sym(fourth_derivative_integrand(E0, E1, E2, E3, E4))
        d4b = sym(lmat * d4s * lmat.T)
        s2 = op_norm_symmetric(d2s)
        b2 = op_norm_symmetric(d2b)
        s4 = op_norm_symmetric(d4s)
        b4 = op_norm_symmetric(d4b)
        e0n = op_norm_rect(E0)
        e1n = op_norm_rect(E1)
        e2n = op_norm_rect(E2)
        max_s2 = max(max_s2, s2)
        max_b2 = max(max_b2, b2)
        max_s4 = max(max_s4, s4)
        max_b4 = max(max_b4, b4)
        max_e = max(max_e, e0n)
        max_e1 = max(max_e1, e1n)
        max_e2 = max(max_e2, e2n)
        rows.append(
            {
                "u": f(u),
                "d2SOp": f(s2),
                "d2BBStarOp": f(b2),
                "d4SOp": f(s4),
                "d4BBStarOp": f(b4),
                "sourceOp": f(e0n),
                "sourceDerivativeOp": f(e1n),
                "sourceSecondDerivativeOp": f(e2n),
            }
        )

    factor = (b - a) * h**2 / 12
    s_error_bound = factor * max_s2
    bgram_error_bound = factor * max_b2
    s_norm = op_norm_symmetric(s_quad)
    bgram_norm = op_norm_symmetric(bgram_quad)
    b_norm = mp.sqrt(bgram_norm)
    s_rel = s_error_bound / s_norm if s_norm else mp.inf
    bgram_rel = bgram_error_bound / bgram_norm if bgram_norm else mp.inf
    # If B B^* has relative error eps and eps<1, singular values change by at
    # most sqrt(eps) in crude relative scale.  The BB* figure is the invariant
    # finite-dimensional quantity; this derived number is intentionally marked
    # as crude.
    crude_b_rel = mp.sqrt(bgram_rel) if bgram_rel >= 0 else mp.inf

    target_s = mp.mpf(args.target_s_rel)
    target_b = mp.mpf(args.target_b_rel)
    print(
        f"Source quadrature error bound basis={args.basis} grid={args.source_grid} "
        f"envelope_grid={args.envelope_grid}",
        flush=True,
    )
    print(f"  ||S_h||={fmt(s_norm, 10)} error_bound={fmt(s_error_bound, 10)} rel={fmt(s_rel, 10)}", flush=True)
    print(
        f"  ||BB*_h||={fmt(bgram_norm, 10)} error_bound={fmt(bgram_error_bound, 10)} "
        f"rel={fmt(bgram_rel, 10)} crude_B_rel={fmt(crude_b_rel, 10)}",
        flush=True,
    )
    print(f"  target S rel={fmt(target_s, 8)} pass={s_rel < target_s}", flush=True)
    print(f"  target crude B rel={fmt(target_b, 8)} pass={crude_b_rel < target_b}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "point": f(mp.mpf(args.point)),
        "orders": args.orders,
        "basis": args.basis,
        "sourceMin": f(a),
        "sourceMax": f(b),
        "sourceGrid": args.source_grid,
        "mesh": f(h),
        "envelopeGrid": args.envelope_grid,
        "trapezoidFactor": f(factor),
        "sourceGramNorm": f(s_norm),
        "sourceGramSecondDerivativeEnvelope": f(max_s2),
        "sourceGramFourthDerivativeEnvelope": f(max_s4),
        "sourceGramErrorBound": f(s_error_bound),
        "sourceGramRelativeErrorBound": f(s_rel),
        "sourceGramTargetRelative": f(target_s),
        "sourceGramPassesTarget": bool(s_rel < target_s),
        "bbStarGramNorm": f(bgram_norm),
        "bOperatorNorm": f(b_norm),
        "bbStarSecondDerivativeEnvelope": f(max_b2),
        "bbStarFourthDerivativeEnvelope": f(max_b4),
        "bbStarErrorBound": f(bgram_error_bound),
        "bbStarRelativeErrorBound": f(bgram_rel),
        "crudeBRelativeErrorFromBBStar": f(crude_b_rel),
        "bTargetRelative": f(target_b),
        "crudeBPassesTarget": bool(crude_b_rel < target_b),
        "maxSourceOp": f(max_e),
        "maxSourceDerivativeOp": f(max_e1),
        "maxSourceSecondDerivativeOp": f(max_e2),
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "envelopeRows": rows,
        "interpretation": (
            "Composite trapezoid error for S=int E^*E and BB*=L S L^*. "
            "The bound is deterministic once the sampled derivative envelope "
            "is replaced by an analytic supremum bound.  S is the finite "
            "Gram representative of T=EE^* on the nonzero source range."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
