#!/usr/bin/env python3
r"""Ball-arithmetic Chebyshev coefficient certificate.

This is the rounded companion to ``source_quadrature_chebyshev_envelope.py``.
It keeps the same finite Green/source evaluator, then wraps each sampled matrix
entry in an explicit interval ball and propagates those balls through the
DCT-I coefficient calculation and the final geometric-tail estimate.

The output certifies the coefficient/tail arithmetic conditional on the stated
sample balls.  A separate Bernstein-ellipse bound supplies the analytic reason
that the geometric tail is legitimate for the explicit Green-source formula.
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
from source_quadrature_chebyshev_envelope import (  # noqa: E402
    cheb_lobatto_nodes,
    matrix_at,
)
from source_quadrature_error_bound import (  # noqa: E402
    op_norm_symmetric,
    trapezoid_source_gram,
)


def iv_const(x):
    x = mp.mpf(x)
    return mp.iv.mpf([x, x])


def iv_ball(x, rel_radius, abs_radius):
    x = mp.mpf(x)
    radius = max(mp.mpf(abs_radius), abs(x) * mp.mpf(rel_radius))
    return mp.iv.mpf([x - radius, x + radius])


def iv_abs_sup(x):
    return max(abs(mp.mpf(x.a)), abs(mp.mpf(x.b)))


def iv_abs_inf(x):
    lo = mp.mpf(x.a)
    hi = mp.mpf(x.b)
    if lo <= 0 <= hi:
        return mp.mpf("0")
    return min(abs(lo), abs(hi))


def iv_cheb_coefficients(values):
    degree = len(values) - 1
    coeffs = []
    two_over_degree = iv_const(mp.mpf("2") / degree)
    half = iv_const(mp.mpf("0.5"))
    one = iv_const(mp.mpf("1"))
    for k in range(degree + 1):
        total = iv_const(0)
        for j, value in enumerate(values):
            weight = half if j == 0 or j == degree else one
            theta = mp.pi * k * j / degree
            cos_theta = mp.iv.cos(mp.iv.mpf([theta, theta]))
            total += weight * value * cos_theta
        coeff = two_over_degree * total
        if k == 0 or k == degree:
            coeff *= half
        coeffs.append(coeff)
    return coeffs


def iv_geometric_tail(coeffs, tail_window):
    degree = len(coeffs) - 1
    start = max(1, degree - tail_window + 1)
    ratios = []
    for k in range(start + 1, degree + 1):
        denominator = iv_abs_inf(coeffs[k - 1])
        if denominator == 0:
            return mp.inf, mp.inf
        ratios.append(iv_abs_sup(coeffs[k]) / denominator)
    q = max(ratios) if ratios else mp.mpf("0")
    if q >= 1:
        return mp.inf, q
    return iv_abs_sup(coeffs[-1]) * q / (1 - q), q


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
    parser.add_argument("--sample-rel-radius", default="1e-55")
    parser.add_argument("--sample-abs-radius", default="1e-60")
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
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--iv-dps", type=int, default=70)
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
    parser.add_argument("--json-out", default="source_quadrature_chebyshev_ball_d32_g129.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace data max_q={max_q} basis={args.basis}", flush=True)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print("Preparing local high block", flush=True)
    base = local_case(args, args.basis, e_derivs)
    print("Preparing endpoint quadrature", flush=True)
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
    sample_rel = mp.mpf(args.sample_rel_radius)
    sample_abs = mp.mpf(args.sample_abs_radius)

    print(
        f"Ball Chebyshev source envelope basis={args.basis} source_grid={args.source_grid} "
        f"degree={args.cheb_degree}",
        flush=True,
    )
    matrices = []
    sampled_op_max = mp.mpf("0")
    max_sample_radius = mp.mpf("0")
    for idx, u in enumerate(nodes):
        mat = matrix_at(args, base, e_derivs, r_nodes, r_weights, u)
        matrices.append(mat)
        sampled_op_max = max(sampled_op_max, op_norm_symmetric(mat))
        for value in mat:
            max_sample_radius = max(max_sample_radius, max(sample_abs, abs(value) * sample_rel))
        if (idx + 1) % max(1, len(nodes) // 8) == 0:
            print(f"  sampled {idx+1:4d}/{len(nodes)} max={fmt(sampled_op_max, 8)}", flush=True)

    dim = matrices[0].rows
    entry_bounds = mp.matrix(dim)
    finite_entry_bounds = mp.matrix(dim)
    tail_entry_bounds = mp.matrix(dim)
    max_tail = mp.mpf("0")
    max_q_tail = mp.mpf("0")
    worst_entry = None
    worst_width = mp.mpf("0")
    tail_failed_entries = []
    for i in range(dim):
        for j in range(dim):
            values = [
                iv_ball(mat[i, j], sample_rel, sample_abs)
                for mat in matrices
            ]
            coeffs = iv_cheb_coefficients(values)
            abs_coeffs_upper = [iv_abs_sup(c) for c in coeffs]
            for c in coeffs:
                worst_width = max(worst_width, mp.mpf(c.b) - mp.mpf(c.a))
            tail, q = iv_geometric_tail(coeffs, args.tail_window)
            if not mp.isfinite(tail):
                tail_failed_entries.append([i, j])
            finite_bound = mp.fsum(abs_coeffs_upper)
            bound = finite_bound + tail
            finite_entry_bounds[i, j] = finite_bound
            tail_entry_bounds[i, j] = tail
            entry_bounds[i, j] = bound
            if tail > max_tail:
                max_tail = tail
                max_q_tail = q
                worst_entry = [i, j]

    row_sum_bound = max(
        mp.fsum(entry_bounds[i, j] for j in range(dim))
        for i in range(dim)
    )
    finite_row_sum_bound = max(
        mp.fsum(finite_entry_bounds[i, j] for j in range(dim))
        for i in range(dim)
    )
    tail_row_sum_bound = max(
        mp.fsum(tail_entry_bounds[i, j] for j in range(dim))
        for i in range(dim)
    )
    frob_bound = mp.sqrt(mp.fsum(entry_bounds[i, j] ** 2 for i in range(dim) for j in range(dim)))
    finite_frob_bound = mp.sqrt(
        mp.fsum(finite_entry_bounds[i, j] ** 2 for i in range(dim) for j in range(dim))
    )
    tail_frob_bound = mp.sqrt(
        mp.fsum(tail_entry_bounds[i, j] ** 2 for i in range(dim) for j in range(dim))
    )
    envelope = min(row_sum_bound, frob_bound)
    finite_envelope = min(finite_row_sum_bound, finite_frob_bound)
    tail_envelope = min(tail_row_sum_bound, tail_frob_bound)
    error_bound = factor * envelope
    rel = error_bound / s_norm if s_norm else mp.inf
    target = mp.mpf(args.target_s_rel)
    pass_target = bool(mp.isfinite(rel) and rel < target)

    print(f"  sampled_op_max={fmt(sampled_op_max, 12)}", flush=True)
    print(f"  ball_envelope={fmt(envelope, 12)} row_sum={fmt(row_sum_bound, 12)}", flush=True)
    print(f"  tail_q_upper={fmt(max_q_tail, 12)} max_tail={fmt(max_tail, 12)}", flush=True)
    print(f"  rel_error={fmt(rel, 12)} target={target} pass={pass_target}", flush=True)

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
        "mpDps": args.dps,
        "ivDps": args.iv_dps,
        "sampleRelativeRadius": f(sample_rel),
        "sampleAbsoluteRadius": f(sample_abs),
        "maxSampleBallRadius": f(max_sample_radius),
        "trapezoidFactor": f(factor),
        "sourceGramNorm": f(s_norm),
        "sampledOperatorEnvelope": f(sampled_op_max),
        "ballChebyshevEnvelope": f(envelope),
        "finiteCoefficientEnvelope": f(finite_envelope),
        "tailEnvelope": f(tail_envelope),
        "rowSumEnvelope": f(row_sum_bound),
        "finiteRowSumEnvelope": f(finite_row_sum_bound),
        "tailRowSumEnvelope": f(tail_row_sum_bound),
        "frobeniusEnvelope": f(frob_bound),
        "finiteFrobeniusEnvelope": f(finite_frob_bound),
        "tailFrobeniusEnvelope": f(tail_frob_bound),
        "maxTailEstimateUpper": f(max_tail),
        "maxTailRatioUpper": f(max_q_tail),
        "worstTailEntry": worst_entry,
        "coefficientMaxWidth": f(worst_width),
        "tailFailedEntries": tail_failed_entries,
        "sourceGramErrorBound": f(error_bound),
        "sourceGramRelativeErrorBound": f(rel),
        "sourceGramTargetRelative": f(target),
        "sourceGramPassesTarget": pass_target,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Ball propagation for the Chebyshev coefficient and final-window "
            "tail arithmetic.  The finite Green-source samples are enclosed "
            "in the displayed sample balls; the Bernstein-ellipse lemma gives "
            "the analytic justification for geometric Chebyshev decay."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
