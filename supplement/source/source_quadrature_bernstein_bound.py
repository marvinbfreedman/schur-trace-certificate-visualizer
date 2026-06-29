#!/usr/bin/env python3
r"""Bernstein-ellipse tail bound for the source Chebyshev certificate.

For ``u`` in the source window [a,b], write ``u=c+r x`` with ``x in [-1,1]``.
If each entry of

    A(u) = d_u^2(E(u)^* E(u))

is analytic in the Bernstein ellipse ``E_rho`` and

    M_rho >= max_i sum_j sup_{z in E_rho} |A_ij(c+r z)|,

then the Chebyshev tail after degree N has row-sum norm at most

    2 M_rho rho^(-(N+1)) / (1 - rho^(-1)).

This script records the geometry and samples the explicit Green-source formula
on the ellipse to produce a numerical M_rho proxy.  The theorem in the notes is
the rigorous part: replacing the sampled M_rho by an interval/complex-ball
supremum gives a fully rigorous Bernstein certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_quadrature_chebyshev_envelope import matrix_at  # noqa: E402


def load_finite_envelope(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("finiteCoefficientEnvelope", "ballChebyshevEnvelope", "chebyshevEnvelope"):
        if key in data:
            return mp.mpf(data[key]), key, data
    raise KeyError(f"{path} has no recognized Chebyshev envelope key")


def bernstein_point(theta, rho):
    w = rho * mp.e ** (1j * theta)
    return (w + 1 / w) / 2


def entry_abs_envelope(mat):
    rows = mat.rows
    cols = mat.cols
    return [[abs(mat[i, j]) for j in range(cols)] for i in range(rows)]


def combine_entry_max(current, mat_abs):
    if current is None:
        return [[value for value in row] for row in mat_abs]
    for i, row in enumerate(mat_abs):
        for j, value in enumerate(row):
            if value > current[i][j]:
                current[i][j] = value
    return current


def matrix_entry_envelope(entry_max):
    row_sum = max(mp.fsum(row) for row in entry_max)
    frob = mp.sqrt(mp.fsum(value * value for row in entry_max for value in row))
    return min(row_sum, frob), row_sum, frob


def formal_pole_rho(a, b, pole):
    center = (a + b) / 2
    radius = (b - a) / 2
    x0 = (pole - center) / radius
    if abs(x0) <= 1:
        return x0, mp.mpf("1")
    root = mp.sqrt(x0 * x0 - 1)
    candidates = [abs(x0 + root), abs(x0 - root)]
    return x0, max(candidates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--finite-envelope-json", default="source_quadrature_chebyshev_ball_d32_g129.json")
    parser.add_argument("--rho", default="2.0")
    parser.add_argument("--ellipse-samples", type=int, default=96)
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
    parser.add_argument("--json-out", default="source_quadrature_bernstein_bound_rho2.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps

    finite_envelope, finite_key, finite_data = load_finite_envelope(args.finite_envelope_json)
    source_norm = mp.mpf(args.source_gram_norm)
    if source_norm <= 0 and "sourceGramNorm" in finite_data:
        source_norm = mp.mpf(finite_data["sourceGramNorm"])
    if source_norm <= 0:
        raise ValueError("source Gram norm must be supplied or present in finite-envelope JSON")

    a = mp.mpf(args.source_min)
    b = mp.mpf(args.source_max)
    center = (a + b) / 2
    radius = (b - a) / 2
    h = (b - a) / (args.source_grid - 1)
    trapezoid_factor = (b - a) * h**2 / 12
    rho = mp.mpf(args.rho)
    tail_factor = 2 * rho ** (-(args.cheb_degree + 1)) / (1 - 1 / rho)
    x0, rho_formula_pole = formal_pole_rho(a, b, mp.mpf("0"))

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    base = local_case(args, args.basis, e_derivs)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    entry_max = None
    max_row_sample = mp.mpf("0")
    print(
        f"Bernstein source envelope basis={args.basis} rho={rho} "
        f"samples={args.ellipse_samples}",
        flush=True,
    )
    for k in range(args.ellipse_samples):
        theta = 2 * mp.pi * k / args.ellipse_samples
        z = bernstein_point(theta, rho)
        u = center + radius * z
        mat = matrix_at(args, base, e_derivs, r_nodes, r_weights, u)
        mat_abs = entry_abs_envelope(mat)
        entry_max = combine_entry_max(entry_max, mat_abs)
        row_sample = max(mp.fsum(row) for row in mat_abs)
        max_row_sample = max(max_row_sample, row_sample)
        if (k + 1) % max(1, args.ellipse_samples // 8) == 0:
            print(f"  ellipse {k+1:4d}/{args.ellipse_samples} row={fmt(max_row_sample, 8)}", flush=True)

    sampled_mrho, sampled_row_mrho, sampled_frob_mrho = matrix_entry_envelope(entry_max)
    bernstein_tail = sampled_mrho * tail_factor
    total_envelope = finite_envelope + bernstein_tail
    error_bound = trapezoid_factor * total_envelope
    rel = error_bound / source_norm
    target = mp.mpf(args.target_s_rel)
    pass_target = bool(rel < target)

    print(f"  finite_envelope={fmt(finite_envelope, 12)} ({finite_key})", flush=True)
    print(f"  sampled_Mrho={fmt(sampled_mrho, 12)} tail_factor={fmt(tail_factor, 12)}", flush=True)
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
        "ellipseSamples": args.ellipse_samples,
        "mappedFormulaPoleAtU0": f(x0),
        "rhoToFormulaPoleU0": f(rho_formula_pole),
        "trapezoidFactor": f(trapezoid_factor),
        "sourceGramNorm": f(source_norm),
        "finiteEnvelopeSource": args.finite_envelope_json,
        "finiteEnvelopeKey": finite_key,
        "finiteCoefficientEnvelope": f(finite_envelope),
        "sampledBernsteinMrho": f(sampled_mrho),
        "sampledBernsteinRowMrho": f(sampled_row_mrho),
        "sampledBernsteinFrobeniusMrho": f(sampled_frob_mrho),
        "bernsteinTailFactor": f(tail_factor),
        "bernsteinTailEnvelope": f(bernstein_tail),
        "bernsteinTotalEnvelope": f(total_envelope),
        "sourceGramErrorBound": f(error_bound),
        "sourceGramRelativeErrorBound": f(rel),
        "sourceGramTargetRelative": f(target),
        "sourceGramPassesTarget": pass_target,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Bernstein-ellipse tail replacement for the Chebyshev source "
            "certificate.  The displayed theorem is rigorous with any true "
            "M_rho bound; this run uses sampled ellipse data as the numerical "
            "M_rho proxy to show the margin size."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
