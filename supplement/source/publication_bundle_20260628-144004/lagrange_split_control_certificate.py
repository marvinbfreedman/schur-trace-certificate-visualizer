#!/usr/bin/env python3
r"""Low/high split for the trace Lagrange residual operator.

The corrected endpoint theorem should split the Lagrange residual control into

  1. a finite low-mode Schur block;
  2. a high-frequency Hardy/commuted-Sturm tail estimate.

For each source point u, this script builds the two residual rows

    ell_0(f)=B_P[h_u,f](s0),
    ell_1(f)=P^*h_u(s0) f(s0),

on the sampled closed trace quotient N=ker R.  If A=K|_N has eigenpairs
``A e_j=lambda_j e_j`` and residual coefficient vectors

    r_j = (ell_0(e_j), ell_1(e_j)) in R^2,

then the exact finite residual control matrix is

    Gamma = sum_j lambda_j^{-1} r_j r_j^T.

The low Schur block is the partial sum j<M.  The high tail is the partial sum
j>=M, and for any moment m>=0 it obeys the finite spectral bound

    Gamma_high(M) <= lambda_M^{-(m+1)} sum_j lambda_j^m r_j r_j^T.

This is the finite analogue of the proposed theorem:

    low Schur block removal + high-frequency Hardy/commuted Sturm trace bound.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_energy_control_certificate import (  # noqa: E402
    boundary_functional_row,
    eval_functional_row,
    make_qargs,
    split_spaces,
)
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    max_eig_or_zero,
    trace_matrix,
)
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_lagrange_adjoint_control import (  # noqa: E402
    adjoint_source_value,
    load_exact_trace,
    trace_green_concomitant_row,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def row_times_mat(row, mat):
    out = mp.matrix(1, mat.cols)
    for j in range(mat.cols):
        out[0, j] = mp.fsum(row[0, i] * mat[i, j] for i in range(mat.rows))
    return out


def residual_coefficients(row_matrix, modes):
    return row_matrix * modes


def two_by_two_from_coeffs(coeffs, avals, start: int, stop: int, inverse_power=1):
    out = mp.matrix(2)
    for j in range(start, stop):
        scale = avals[j] ** (-mp.mpf(inverse_power))
        for a in range(2):
            for b in range(2):
                out[a, b] += scale * coeffs[a, j] * coeffs[b, j]
    return (out + out.T) / 2


def moment_matrix(coeffs, avals, moment: float):
    out = mp.matrix(2)
    p = mp.mpf(moment)
    for j, lam in enumerate(avals):
        scale = lam**p
        for a in range(2):
            for b in range(2):
                out[a, b] += scale * coeffs[a, j] * coeffs[b, j]
    return (out + out.T) / 2


def matrix_entries(mat):
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def compute_case(args, t, polys, N, modes, avals, e_derivs, r_nodes, r_weights):
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h_derivs = source_derivatives(
        c,
        s0,
        t,
        args.jet_order,
        r_nodes,
        r_weights,
    )
    pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
    brow = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
    brow_poly = boundary_functional_row(polys, s0, brow)
    eval_poly = eval_functional_row(polys, s0, pstar)

    residual_rows = mp.matrix(2, len(polys))
    for j in range(len(polys)):
        residual_rows[0, j] = brow_poly[0, j]
        residual_rows[1, j] = eval_poly[0, j]
    coeffs = residual_coefficients(residual_rows, N * modes)

    full = two_by_two_from_coeffs(coeffs, avals, 0, len(avals), inverse_power=1)
    full_const = max_eig_or_zero(full)
    splits = []
    for cutoff in parse_ints(args.cutoffs):
        cutoff = min(max(0, cutoff), len(avals))
        low = two_by_two_from_coeffs(coeffs, avals, 0, cutoff, inverse_power=1)
        high = two_by_two_from_coeffs(coeffs, avals, cutoff, len(avals), inverse_power=1)
        lambda_cut = avals[cutoff] if cutoff < len(avals) else mp.inf
        moment_bounds = []
        for moment in parse_floats(args.moments):
            mm = moment_matrix(coeffs, avals, moment)
            if cutoff < len(avals):
                bound = (lambda_cut ** (-(mp.mpf(moment) + 1))) * max_eig_or_zero(mm)
            else:
                bound = mp.mpf("0")
            moment_bounds.append(
                {
                    "moment": moment,
                    "momentOp": f(max_eig_or_zero(mm)),
                    "highBound": f(bound),
                    "highBoundFrac": f(bound / full_const if full_const else mp.mpf("0")),
                }
            )
        splits.append(
            {
                "cutoff": cutoff,
                "lambdaCut": f(lambda_cut) if cutoff < len(avals) else None,
                "lowConst": f(max_eig_or_zero(low)),
                "highConst": f(max_eig_or_zero(high)),
                "lowFrac": f(max_eig_or_zero(low) / full_const if full_const else mp.mpf("0")),
                "highFrac": f(max_eig_or_zero(high) / full_const if full_const else mp.mpf("0")),
                "lowMatrix": matrix_entries(low),
                "highMatrix": matrix_entries(high),
                "momentBounds": moment_bounds,
            }
        )

    mode_rows = []
    for j, lam in enumerate(avals):
        rr = mp.matrix(2, 1)
        rr[0] = coeffs[0, j]
        rr[1] = coeffs[1, j]
        single = (rr * rr.T) / lam
        mode_rows.append(
            {
                "mode": j,
                "lambda": f(lam),
                "singleConst": f(max_eig_or_zero(single)),
                "singleFrac": f(max_eig_or_zero(single) / full_const if full_const else mp.mpf("0")),
                "boundaryCoeff": f(coeffs[0, j]),
                "adjointCoeff": f(coeffs[1, j]),
            }
        )

    return {
        "t": f(t),
        "pStarH": f(pstar),
        "fullConst": f(full_const),
        "fullNorm": f(mp.sqrt(max(mp.mpf("0"), full_const))),
        "fullMatrix": matrix_entries(full),
        "splits": splits,
        "modes": mode_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=16)
    parser.add_argument("--quad", type=int, default=22)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=10)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoffs", default="1 2 3 4 5")
    parser.add_argument("--moments", default="0 1 2 3")
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
    parser.add_argument("--json-out", default="lagrange_split_control_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    qargs = make_qargs(args)
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    _centers, R = trace_matrix(polys, qargs)
    N, _U, rank, nullity = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    modes = columns(avecs_all, keep)

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]

    print(
        f"Lagrange split control basis={args.basis} constraints={args.constraints} "
        f"nullity={nullity} positive={len(avals)}",
        flush=True,
    )
    print(
        f"  lambda range=[{fmt(avals[0], 8)}, {fmt(avals[-1], 8)}]",
        flush=True,
    )
    print("  t       full_norm   high@2    high@3    high@4    high@5", flush=True)

    rows = []
    for t in t_values:
        row = compute_case(args, t, polys, N, modes, avals, e_derivs, r_nodes, r_weights)
        rows.append(row)
        by_cut = {item["cutoff"]: item for item in row["splits"]}
        print(
            f"  {fmt(t, 6):>6} {fmt(mp.mpf(row['fullNorm']), 8):>11} "
            f"{fmt(mp.mpf(by_cut.get(2, {'highFrac': 0})['highFrac']), 8):>8} "
            f"{fmt(mp.mpf(by_cut.get(3, {'highFrac': 0})['highFrac']), 8):>8} "
            f"{fmt(mp.mpf(by_cut.get(4, {'highFrac': 0})['highFrac']), 8):>8} "
            f"{fmt(mp.mpf(by_cut.get(5, {'highFrac': 0})['highFrac']), 8):>8}",
            flush=True,
        )

    data = {
        "s0": f(mp.mpf(args.s0)),
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "lambdaMin": f(avals[0]),
        "lambdaMax": f(avals[-1]),
        "cutoffs": parse_ints(args.cutoffs),
        "moments": parse_floats(args.moments),
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
