#!/usr/bin/env python3
r"""Boundary-row Green representer certificate.

The remaining closed-trace Hardy/Green row is

    B_P[h_u,f](s0) = sum_{j=0}^7 b_j(u) f^(j)(s0).

Thus, as with the adjoint-evaluation row, the continuum representer is a
finite combination of fixed objects:

    g_u^bdry = sum_j b_j(u) k_j^{hi},

where k_j^{hi} is the A-Green representer of the jet evaluation
f -> f^(j)(s0) on the closed-trace high block.  The differentiated boundary
row is controlled by replacing b_j(u) by b'_j(u).

This script constructs the fixed jet-representer Gram matrix

    G_ij = <k_i^{hi}, k_j^{hi}>_A,

verifies the factorization of the boundary row, and gives compact source-window
envelopes for ||b(u)||_G and ||b'(u)||_G using b''.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_eval_representer_certificate import (  # noqa: E402
    energy_representer,
    full_coeffs,
    source_derivatives_du_order,
)
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    poly_derivative_value,
    trace_matrix,
)
from trace_lagrange_adjoint_control import (  # noqa: E402
    load_exact_trace,
    trace_green_concomitant_row,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def jet_eval_row(polys, s0, deriv):
    row = mp.matrix(1, len(polys))
    for i, poly in enumerate(polys):
        row[0, i] = poly_derivative_value(poly, s0, deriv)
    return row


def row_from_coeffs(polys, s0, coeffs):
    row = mp.matrix(1, len(polys))
    for deriv, coeff in enumerate(coeffs):
        for i, poly in enumerate(polys):
            row[0, i] += coeff * poly_derivative_value(poly, s0, deriv)
    return row


def gram_quadratic(G, coeffs):
    return mp.fsum(
        coeffs[i] * G[i, j] * coeffs[j]
        for i in range(len(coeffs))
        for j in range(len(coeffs))
    )


def gram_norm(G, coeffs):
    return mp.sqrt(max(mp.mpf("0"), gram_quadratic(G, coeffs)))


def matrix_to_lists(mat):
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


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
    parser.add_argument("--json-out", default="boundary_row_representer_certificate.json")
    args = parser.parse_args()

    if args.jet_order < 2:
        raise SystemExit("--jet-order must be at least 2")
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

    order = args.jet_order - 1
    jet_records = []
    jet_reps = []
    for deriv in range(order):
        row_n = jet_eval_row(polys, s0, deriv) * N
        rep_n, norm2, rel = energy_representer(row_n, A, high_a_modes, high_avals)
        jet_reps.append(rep_n)
        coeffs = full_coeffs(N, rep_n)
        jet_records.append(
            {
                "deriv": deriv,
                "norm2": f(norm2),
                "norm": f(mp.sqrt(max(mp.mpf("0"), norm2))),
                "rangeRelativeDefect": f(rel),
                "coeffMaxAbs": f(max(abs(c) for c in coeffs) if coeffs else mp.mpf("0")),
            }
        )

    G = mp.matrix(order)
    for i in range(order):
        for j in range(order):
            G[i, j] = (jet_reps[i].T * A * jet_reps[j])[0, 0]

    vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    source_min = mp.mpf(args.source_min)
    source_max = mp.mpf(args.source_max)
    step = (source_max - source_min) / (args.source_grid - 1)
    radius = step / 2
    c = mp.pi

    rows = []
    max_b_norm = mp.mpf("0")
    max_db_norm = mp.mpf("0")
    max_d2b_norm = mp.mpf("0")
    max_factor_defect = mp.mpf("0")
    for i in range(args.source_grid):
        u = source_min + i * step
        h0 = source_derivatives_du_order(c, s0, u, args.jet_order, 0, r_nodes, r_weights)
        h1 = source_derivatives_du_order(c, s0, u, args.jet_order, 1, r_nodes, r_weights)
        h2 = source_derivatives_du_order(c, s0, u, args.jet_order, 2, r_nodes, r_weights)
        b = trace_green_concomitant_row(e_derivs, h0, args.jet_order)
        db = trace_green_concomitant_row(e_derivs, h1, args.jet_order)
        d2b = trace_green_concomitant_row(e_derivs, h2, args.jet_order)

        predicted = gram_quadratic(G, b)
        predicted_d = gram_quadratic(G, db)
        row_n = row_from_coeffs(polys, s0, b) * N
        drow_n = row_from_coeffs(polys, s0, db) * N
        _rep, actual, _rel = energy_representer(row_n, A, high_a_modes, high_avals)
        _drep, actual_d, _drel = energy_representer(drow_n, A, high_a_modes, high_avals)
        scale = max(mp.mpf("1"), abs(predicted), abs(actual), abs(predicted_d), abs(actual_d))
        factor_defect = max(abs(actual - predicted), abs(actual_d - predicted_d)) / scale
        max_factor_defect = max(max_factor_defect, factor_defect)

        b_norm = gram_norm(G, b)
        db_norm = gram_norm(G, db)
        d2b_norm = gram_norm(G, d2b)
        max_b_norm = max(max_b_norm, b_norm)
        max_db_norm = max(max_db_norm, db_norm)
        max_d2b_norm = max(max_d2b_norm, d2b_norm)
        rows.append(
            {
                "u": f(u),
                "boundaryNorm2": f(actual),
                "dBoundaryNorm2": f(actual_d),
                "bNorm": f(b_norm),
                "dBNorm": f(db_norm),
                "d2BNorm": f(d2b_norm),
                "factorRelativeDefect": f(factor_defect),
                "coefficients": [f(x) for x in b],
                "dCoefficients": [f(x) for x in db],
                "d2Coefficients": [f(x) for x in d2b],
            }
        )

    b_cover = max_b_norm + radius * max_db_norm
    db_cover = max_db_norm + radius * max_d2b_norm
    boundary_cover = b_cover * b_cover
    dboundary_cover = db_cover * db_cover

    print(
        f"Boundary-row representer certificate basis={args.basis} "
        f"constraints={args.constraints} cutoff={cutoff}",
        flush=True,
    )
    print(f"  max jet range defect    = {fmt(max(mp.mpf(r['rangeRelativeDefect']) for r in jet_records), 8)}", flush=True)
    print(f"  max factor defect       = {fmt(max_factor_defect, 8)}", flush=True)
    print(
        f"  cover ||b||_G={fmt(b_cover, 10)} ||b'||_G={fmt(db_cover, 10)} "
        f"bdry={fmt(boundary_cover, 10)} dB={fmt(dboundary_cover, 10)}",
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
        "traceLambda0": f(vals[0]),
        "traceLambda1": f(vals[1]),
        "sourceMin": f(source_min),
        "sourceMax": f(source_max),
        "sourceGrid": args.source_grid,
        "mesh": f(step),
        "meshRadius": f(radius),
        "maxJetRangeRelativeDefect": f(max(mp.mpf(r["rangeRelativeDefect"]) for r in jet_records)),
        "maxFactorRelativeDefect": f(max_factor_defect),
        "maxBNormGrid": f(max_b_norm),
        "maxDBNormGrid": f(max_db_norm),
        "maxD2BNormGrid": f(max_d2b_norm),
        "coverBNorm": f(b_cover),
        "coverDBNorm": f(db_cover),
        "coverBoundaryNorm2": f(boundary_cover),
        "coverDBoundaryNorm2": f(dboundary_cover),
        "jetRepresenters": jet_records,
        "jetGram": matrix_to_lists(G),
        "rows": rows,
        "theoremTemplate": (
            "The boundary-row continuum representer is a finite sum of fixed "
            "closed-trace high-block jet representers. Uniform control follows "
            "from boundedness of the fixed jet Gram matrix and compact-window "
            "bounds for the coefficient row b(u) and its derivative."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
