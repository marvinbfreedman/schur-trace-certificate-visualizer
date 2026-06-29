#!/usr/bin/env python3
r"""Chebyshev-tail envelope for d_u^2(E(u)^*E(u)).

Raw interval arithmetic over the Green-source formula is far too loose because
of dependency in the exponential/Taylor representation.  This script uses the
analytic Chebyshev route instead:

1. sample the matrix entries of A(u)=d_u^2(E(u)^*E(u)) at Chebyshev-Lobatto
   nodes on the source window;
2. compute Chebyshev coefficients entrywise;
3. bound sup_u ||A(u)|| by a row-sum bound from coefficient l1 norms plus a
   geometric tail estimate from the final coefficients.

The output is a formal certificate modulo the displayed Chebyshev tail
hypothesis.  To turn it into a fully machine-rigorous certificate, replace the
floating coefficient/tail arithmetic by ball arithmetic or prove the same tail
bound from a Bernstein ellipse estimate.
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
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_quadrature_error_bound import (  # noqa: E402
    op_norm_symmetric,
    source_row,
    trapezoid_source_gram,
)
from trace_active_derivative_rank import derivative_matrix  # noqa: E402


def cheb_lobatto_nodes(a, b, degree):
    center = (a + b) / 2
    radius = (b - a) / 2
    return [center + radius * mp.cos(mp.pi * j / degree) for j in range(degree + 1)]


def jf(x):
    """JSON-safe float conversion."""
    y = float(x)
    if not mp.isfinite(x):
        return None
    return y


def cheb_coefficients(values):
    """DCT-I coefficients for the Chebyshev interpolant."""
    degree = len(values) - 1
    coeffs = []
    for k in range(degree + 1):
        total = mp.mpf("0")
        for j, value in enumerate(values):
            weight = mp.mpf("0.5") if j == 0 or j == degree else mp.mpf("1")
            total += weight * value * mp.cos(mp.pi * k * j / degree)
        coeff = (2 / degree) * total
        if k == 0 or k == degree:
            coeff /= 2
        coeffs.append(coeff)
    return coeffs


def geometric_tail(abs_coeffs, tail_window):
    degree = len(abs_coeffs) - 1
    start = max(1, degree - tail_window + 1)
    ratios = []
    for k in range(start + 1, degree + 1):
        prev = abs_coeffs[k - 1]
        if prev:
            ratios.append(abs_coeffs[k] / prev)
    q = max(ratios) if ratios else mp.mpf("0")
    if q >= 1:
        return mp.inf, q
    return abs_coeffs[-1] * q / (1 - q), q


def cross(A, B):
    return A.T * B


def d2_source_gram(E0, E1, E2):
    return sym(cross(E2, E0) + 2 * cross(E1, E1) + cross(E0, E2))


def matrix_at(args, base, e_derivs, r_nodes, r_weights, u):
    E0 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 0)
    E1 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 1)
    E2 = source_row(args, base, e_derivs, r_nodes, r_weights, u, 2)
    return d2_source_gram(E0, E1, E2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=129)
    parser.add_argument("--source-gram-norm", default="0")
    parser.add_argument("--cheb-degree", type=int, default=32)
    parser.add_argument("--tail-window", type=int, default=8)
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
    parser.add_argument("--json-out", default="source_quadrature_chebyshev_envelope.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    if mp.mpf(args.source_gram_norm) > 0:
        s_norm = mp.mpf(args.source_gram_norm)
    else:
        qargs = SimpleNamespace(**vars(args))
        qargs.source_grid = args.source_grid
        _nodes, _weights, s_quad = trapezoid_source_gram(qargs, base, e_derivs, args.source_grid)
        s_norm = op_norm_symmetric(s_quad)

    a = mp.mpf(args.source_min)
    b = mp.mpf(args.source_max)
    h = (b - a) / (args.source_grid - 1)
    factor = (b - a) * h**2 / 12
    nodes = cheb_lobatto_nodes(a, b, args.cheb_degree)

    print(
        f"Chebyshev source envelope basis={args.basis} source_grid={args.source_grid} "
        f"degree={args.cheb_degree}",
        flush=True,
    )
    matrices = []
    sampled_op_max = mp.mpf("0")
    for idx, u in enumerate(nodes):
        mat = matrix_at(args, base, e_derivs, r_nodes, r_weights, u)
        matrices.append(mat)
        sampled_op_max = max(sampled_op_max, op_norm_symmetric(mat))
        if (idx + 1) % max(1, len(nodes) // 8) == 0:
            print(f"  sampled {idx+1:4d}/{len(nodes)} max={fmt(sampled_op_max, 8)}", flush=True)

    dim = matrices[0].rows
    entry_bounds = mp.matrix(dim)
    max_tail = mp.mpf("0")
    max_q = mp.mpf("0")
    worst_entry = None
    for i in range(dim):
        for j in range(dim):
            values = [mat[i, j] for mat in matrices]
            coeffs = cheb_coefficients(values)
            abs_coeffs = [abs(c) for c in coeffs]
            tail, q = geometric_tail(abs_coeffs, args.tail_window)
            bound = mp.fsum(abs_coeffs) + tail
            entry_bounds[i, j] = bound
            if tail > max_tail:
                max_tail = tail
                max_q = q
                worst_entry = [i, j]

    row_sum_bound = max(
        mp.fsum(entry_bounds[i, j] for j in range(dim))
        for i in range(dim)
    )
    frob_bound = mp.sqrt(mp.fsum(entry_bounds[i, j] ** 2 for i in range(dim) for j in range(dim)))
    envelope = min(row_sum_bound, frob_bound)
    error_bound = factor * envelope
    rel = error_bound / s_norm if s_norm else mp.inf
    target = mp.mpf(args.target_s_rel)
    print(f"  sampled_op_max={fmt(sampled_op_max, 12)}", flush=True)
    print(f"  cheb_envelope={fmt(envelope, 12)} row_sum={fmt(row_sum_bound, 12)}", flush=True)
    print(f"  rel_error={fmt(rel, 12)} target={target} pass={rel < target}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "sourceMin": f(a),
        "sourceMax": f(b),
        "sourceGrid": args.source_grid,
        "sourceGramNormProvided": f(mp.mpf(args.source_gram_norm)),
        "chebDegree": args.cheb_degree,
        "tailWindow": args.tail_window,
        "trapezoidFactor": f(factor),
        "sourceGramNorm": f(s_norm),
        "sampledOperatorEnvelope": f(sampled_op_max),
        "chebyshevEnvelope": jf(envelope),
        "rowSumEnvelope": jf(row_sum_bound),
        "frobeniusEnvelope": jf(frob_bound),
        "maxTailEstimate": jf(max_tail),
        "maxTailRatio": jf(max_q),
        "worstTailEntry": worst_entry,
        "sourceGramErrorBound": jf(error_bound),
        "sourceGramRelativeErrorBound": jf(rel),
        "sourceGramTargetRelative": f(target),
        "sourceGramPassesTarget": bool(rel < target),
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Chebyshev coefficient envelope for d_u^2(E^*E).  The pass is "
            "conditional on the displayed geometric tail estimate; making it "
            "fully rigorous requires ball arithmetic or a Bernstein ellipse "
            "bound for the same coefficients."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
