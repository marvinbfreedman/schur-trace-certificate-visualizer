#!/usr/bin/env python3
r"""Test whether source integration-by-parts boundary rows lie in the trace ideal.

For a scalar source h(s), Q integrations by parts in s produce the boundary row

    B_Q[h](f) = sum_{j=0}^{Q-1} (-1)^j h^(Q-1-j)(s) f^(j)(s)

up to an overall sign convention.  The moving closed trace kills this boundary
row if it lies in the span of

    D_s^q Lambda_s,    0 <= q <= Q-9,

because Lambda_s has order 8.  This script computes h(s)=K_B(s,t) for the
endpoint B-model using formal Taylor coefficients in s, and tests the row-space
membership against the sampled differential closure of Lambda_s.

This is a diagnostic for the proposed s-boundary concomitant closure.  A small
residual would close the algebraic membership step for this candidate
commutation; a large residual says the commuted source needs a different
operator/weight than raw repeated integration by parts.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at, dot  # noqa: E402
from endpoint_kb_confluent_mp import endpoint_series  # noqa: E402
from lambda_differential_closure import (  # noqa: E402
    coefficient_derivatives,
    finite_diff_weights,
    trace_derivative_row,
)
from quotient_factorization_mp import endpoint_b_quadrature, endpoint_b_vw  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def signed_defect(s, args):
    vals, neg, vec = defect_at(
        s,
        args.jet_order,
        mp.mpf(args.defect_rmax),
        args.defect_order,
        mp.mpf(args.tol),
    )
    if vec[0] > 0:
        vec = [-x for x in vec]
    return vals, neg, vec


def sample_trace_field(args):
    s_min = mp.mpf(args.s_min)
    s_max = mp.mpf(args.s_max)
    h = (s_max - s_min) / (args.samples - 1)
    centers = [s_min + i * h for i in range(args.samples)]
    vecs = []
    prev = None
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    for s in centers:
        vals, neg, vec = signed_defect(s, args)
        if prev is not None and dot(prev, vec) < 0:
            vec = [-x for x in vec]
        if prev is not None:
            min_cos = min(min_cos, abs(dot(prev, vec)))
        prev = vec
        vecs.append(vec)
        min_gap = min(min_gap, vals[1] - vals[0])
    return centers, vecs, min_gap, min_cos


def closest_center_index(centers, s0):
    return min(range(len(centers)), key=lambda idx: abs(centers[idx] - s0))


def source_taylor_coeffs(c, s0, t, q_order, r_nodes, r_weights):
    """Return Taylor coefficients h_i with h(s0+z)=sum h_i z^i."""
    coeffs = [mp.mpf("0") for _ in range(q_order)]
    for r, weight in zip(r_nodes, r_weights):
        x = mp.e**r
        tau = x - 1
        _, _, _, vb_series, wb_series = endpoint_series(c, s0, r, q_order)
        vt, wt = endpoint_b_vw(t, r, c)
        layer_weight = weight * x ** mp.mpf("2.5")
        for i in range(q_order):
            coeffs[i] += layer_weight * (
                r * vb_series[i] * vt
                + mp.mpf("0.5") * (vb_series[i] * wt + wb_series[i] * vt)
            )
    return coeffs


def source_derivatives(c, s0, t, max_deriv, r_nodes, r_weights):
    coeffs = source_taylor_coeffs(c, s0, t, max_deriv + 1, r_nodes, r_weights)
    return [mp.factorial(i) * coeffs[i] for i in range(max_deriv + 1)]


def boundary_row_from_derivatives(derivs, q_integrations):
    row = []
    for j in range(q_integrations):
        row.append(((-1) ** j) * derivs[q_integrations - 1 - j])
    return row


def pad_row(row, cols):
    return row + [mp.mpf("0") for _ in range(cols - len(row))]


def solve_membership(target, generators):
    cols = len(target)
    if not generators:
        norm = mp.sqrt(mp.fsum(x * x for x in target))
        return [], norm, norm, mp.mpf("1")
    mat = mp.matrix(cols, len(generators))
    for j, row in enumerate(generators):
        prow = pad_row(row, cols)
        for i, value in enumerate(prow):
            mat[i, j] = value
    b = mp.matrix(cols, 1)
    for i, value in enumerate(target):
        b[i] = value
    gram = mat.T * mat
    rhs = mat.T * b
    try:
        alpha = mp.lu_solve(gram, rhs)
    except ZeroDivisionError:
        alpha = mp.qr_solve(mat, b)[0]
    residual = b - mat * alpha
    res_norm = mp.sqrt(mp.fsum(residual[i] * residual[i] for i in range(cols)))
    target_norm = mp.sqrt(mp.fsum(x * x for x in target))
    rel = res_norm / target_norm if target_norm else mp.mpf("0")
    return [alpha[i] for i in range(alpha.rows)], res_norm, target_norm, rel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s-min", default="0.02")
    parser.add_argument("--s-max", default="0.545")
    parser.add_argument("--samples", type=int, default=17)
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--q", default="9 10 11 12 13")
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--stencil-width", type=int, default=9)
    parser.add_argument("--defect-order", type=int, default=55)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="source_concomitant_membership.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    q_values = [int(piece) for piece in args.q.replace(",", " ").split()]
    max_q = max(q_values)
    max_lambda_q = max_q - args.jet_order
    if max_lambda_q < 0:
        raise SystemExit("all q values must be at least jet_order")

    centers, vecs, min_gap, min_cos = sample_trace_field(args)
    idx = closest_center_index(centers, s0)
    e_derivs = coefficient_derivatives(
        centers, vecs, idx, max_lambda_q, args.stencil_width
    )
    actual_s0 = centers[idx]

    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    rows = []

    print(
        f"Source concomitant membership s0={s0} nearest={actual_s0} "
        f"q={q_values}",
        flush=True,
    )
    print(f"  trace min_gap={fmt(min_gap, 12)} min_cos={fmt(min_cos, 12)}", flush=True)
    print("  t        Q  generators  rel_residual  alpha_norm  target_norm", flush=True)

    for t in t_values:
        derivs = source_derivatives(c, actual_s0, t, max_q - 1, r_nodes, r_weights)
        for q_int in q_values:
            target = boundary_row_from_derivatives(derivs, q_int)
            generators = [
                trace_derivative_row(e_derivs, q, args.jet_order)
                for q in range(q_int - args.jet_order + 1)
            ]
            alpha, res_norm, target_norm, rel = solve_membership(target, generators)
            alpha_norm = mp.sqrt(mp.fsum(a * a for a in alpha)) if alpha else mp.mpf("0")
            rows.append(
                {
                    "t": f(t),
                    "qIntegrations": q_int,
                    "generators": len(generators),
                    "relativeResidual": f(rel),
                    "residualNorm": f(res_norm),
                    "targetNorm": f(target_norm),
                    "alphaNorm": f(alpha_norm),
                    "alphaMaxAbs": f(max(abs(a) for a in alpha)) if alpha else 0.0,
                }
            )
            print(
                f"  {fmt(t, 6):>8} {q_int:2d} {len(generators):10d} "
                f"{fmt(rel, 8):>12} {fmt(alpha_norm, 8):>11} {fmt(target_norm, 8):>11}",
                flush=True,
            )

    data = {
        "sMin": f(mp.mpf(args.s_min)),
        "sMax": f(mp.mpf(args.s_max)),
        "requestedS0": f(s0),
        "actualS0": f(actual_s0),
        "samples": args.samples,
        "jetOrder": args.jet_order,
        "stencilWidth": args.stencil_width,
        "qValues": q_values,
        "tValues": [f(t) for t in t_values],
        "traceMinGap": f(min_gap),
        "traceMinConsecutiveCos": f(min_cos),
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
