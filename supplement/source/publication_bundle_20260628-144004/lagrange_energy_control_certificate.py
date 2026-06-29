#!/usr/bin/env python3
r"""Finite energy control for the trace Lagrange residuals.

The local Lagrange identity leaves two objects to control on the closed trace
space:

    B_P[h_u,f](s0),             f(s0) P^*h_u(s0).

This script measures them as finite Galerkin functionals on

    N = ker R,    Rf=(Lambda_a f)_a,

with the positive endpoint energy

    A = K|_N

as the control norm.  For a row ell on the polynomial basis, the sharp finite
constant is

    ||ell||_{A^{-1}}^2 = ell_N A^+ ell_N^T.

This is not a continuum proof; it is the finite certificate for the exact
control theorem that now has to be proved analytically by a Hardy/Sturm trace
estimate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import (  # noqa: E402
    columns,
    gram_matrix,
    poly_derivative_value,
    trace_matrix,
)
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_lagrange_adjoint_control import (  # noqa: E402
    adjoint_source_value,
    load_exact_trace,
    trace_green_concomitant_row,
)
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def make_qargs(args):
    return SimpleNamespace(
        model="endpoint_b",
        kind="tilde3",
        omega=args.omega,
        L=args.L,
        basis=args.basis,
        quad=args.quad,
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=args.constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
        margin=args.margin,
        dps=args.dps,
    )


def split_spaces(R, rank_tol_text: str):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx), len(u_idx), len(n_idx)


def positive_inverse(mat, tol_text: str):
    vals, vecs = mp.eigsy((mat + mat.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > mp.mpf(tol_text)]
    inv = mp.matrix(mat.rows)
    if keep:
        vkeep = columns(vecs, keep)
        diag = mp.matrix(len(keep))
        for j, idx in enumerate(keep):
            diag[j, j] = 1 / vals[idx]
        inv = vkeep * diag * vkeep.T
    return vals, keep, inv


def row_times_mat(row, mat):
    out = mp.matrix(1, mat.cols)
    for j in range(mat.cols):
        out[0, j] = mp.fsum(row[i] * mat[i, j] for i in range(mat.rows))
    return out


def row_energy_norm2(row, N, Aplus):
    row_n = row_times_mat(row, N)
    return (row_n * Aplus * row_n.T)[0, 0]


def boundary_functional_row(polys, s0, brow):
    row = mp.matrix(1, len(polys))
    for i, poly in enumerate(polys):
        row[0, i] = mp.fsum(
            brow[j] * poly_derivative_value(poly, s0, j)
            for j in range(len(brow))
        )
    return row


def eval_functional_row(polys, s0, scale):
    row = mp.matrix(1, len(polys))
    for i, poly in enumerate(polys):
        row[0, i] = scale * poly_derivative_value(poly, s0, 0)
    return row


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
    parser.add_argument("--json-out", default="lagrange_energy_control_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]

    qargs = make_qargs(args)
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    _centers, R = trace_matrix(polys, qargs)
    N, U, rank, nullity = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    avals, keep, Aplus = positive_inverse(A, args.psd_tol)

    vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    rows = []

    print(
        f"Lagrange energy control basis={args.basis} constraints={args.constraints} "
        f"nullity={nullity} s0={s0}",
        flush=True,
    )
    print(
        f"  A lambda range=[{fmt(avals[0], 8)}, {fmt(avals[-1], 8)}] "
        f"positive={len(keep)} trace_rank={rank}",
        flush=True,
    )
    print("  t       ||B||_A*     ||eval P*h||_A*  combined", flush=True)

    for t in t_values:
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
        b2 = row_energy_norm2(brow_poly, N, Aplus)
        e2 = row_energy_norm2(eval_poly, N, Aplus)
        combined = mp.sqrt(max(mp.mpf("0"), b2)) + mp.sqrt(max(mp.mpf("0"), e2))
        rows.append(
            {
                "t": f(t),
                "pStarH": f(pstar),
                "boundaryNorm2": f(b2),
                "boundaryNorm": f(mp.sqrt(max(mp.mpf("0"), b2))),
                "adjointEvalNorm2": f(e2),
                "adjointEvalNorm": f(mp.sqrt(max(mp.mpf("0"), e2))),
                "combinedNorm": f(combined),
            }
        )
        print(
            f"  {fmt(t, 6):>6} {fmt(mp.sqrt(max(mp.mpf('0'), b2)), 8):>12} "
            f"{fmt(mp.sqrt(max(mp.mpf('0'), e2)), 8):>16} {fmt(combined, 8):>10}",
            flush=True,
        )

    data = {
        "s0": f(s0),
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(keep),
        "lambdaMin": f(avals[0]),
        "lambdaMax": f(avals[-1]),
        "traceLambda0": f(vals[0]),
        "traceLambda1": f(vals[1]),
        "rows": rows,
        "interpretation": (
            "Finite A^{-1} norms are the Galerkin constants for controlling "
            "the endpoint concomitant and local adjoint source by endpoint energy."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
