#!/usr/bin/env python3
r"""Auxiliary high-block regularizer certificate.

Two possible regularizers are compared on H_M cap ker R:

1. Sobolev graph regularizer

       W_10(f)=sum_{r=0}^{10} int |f^(r)|^2.

   It controls local residual rows by Sobolev trace, but may be too large to
   absorb into the positive block A.

2. Source-range residual regularizer

       S_src(f)=sum_t ||E_t f||^2,

   where E_t is the two-row Lagrange residual
   (B_P[h_t,f](s0), P^*h_t(s0)f(s0)).

   This is the minimal finite auxiliary object for the high-block Schur tail:
   it controls the residual by definition, and absorbability is exactly

       S_src(f) <= eta <A f,f>.

The script reports eta for both routes on the same Galerkin/closed-trace block.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_commuted_kernel_energy import positive_inverse_constant  # noqa: E402
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import (  # noqa: E402
    graph_constant,
    residual_rows_for,
    sobolev_matrix,
)
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    max_eig_or_zero,
    trace_matrix,
)
from trace_lagrange_adjoint_control import load_exact_trace  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_floats(text: str) -> list[mp.mpf]:
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def stack_residual_rows(args, polys):
    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    blocks = []
    for t in parse_floats(args.t_values):
        rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, t)
        blocks.append((t, pstar, rows))
    out = mp.matrix(2 * len(blocks), len(polys))
    for block_index, (_t, _pstar, rows) in enumerate(blocks):
        for i in range(2):
            for j in range(len(polys)):
                out[2 * block_index + i, j] = rows[i, j]
    return blocks, out


def generalized_row_constant(row_matrix, modes, denom, tol):
    if modes.cols == 0:
        return mp.mpf("0")
    coeffs = row_matrix * modes
    restricted = modes.T * denom * modes
    _vals, _vecs = mp.eigsy((restricted + restricted.T) / 2, eigvals_only=False)
    _const, _denom_vals = positive_inverse_constant(
        coeffs.T * coeffs,
        restricted,
        tol,
    )
    return _const


def direct_absorption_constant(row_matrix, high_a_modes, avals):
    if high_a_modes.cols == 0:
        return mp.mpf("0")
    # high_a_modes are A-eigenvectors in N-coordinates.  A is diagonal with
    # entries avals on this basis.
    coeffs = row_matrix * high_a_modes
    control = coeffs * coeffs.T
    # Nonzero eigenvalues of C diag(1/lambda) C^T equal those of the
    # A^{-1/2} source Gram on the high block.
    scaled = mp.matrix(coeffs.rows)
    for i in range(coeffs.rows):
        for j in range(coeffs.rows):
            scaled[i, j] = mp.fsum(
                coeffs[i, k] * coeffs[j, k] / avals[k]
                for k in range(len(avals))
            )
    return max_eig_or_zero((scaled + scaled.T) / 2), max_eig_or_zero(control)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
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
    parser.add_argument("--sobolev-order", type=int, default=10)
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
    parser.add_argument("--tol", default="1e-40")
    parser.add_argument("--json-out", default="aux_regularizer_certificate.json")
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
    a_modes = columns(avecs_all, keep)
    full_modes = N * a_modes

    cutoff = min(args.cutoff, len(avals))
    high_modes = columns(full_modes, list(range(cutoff, len(avals))))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]

    blocks, residual_rows = stack_residual_rows(args, polys)
    residual_on_n = residual_rows * N
    wmat_full = sobolev_matrix(polys, mp.mpf(args.L), args.sobolev_order, top_only=False)
    wmat_n = N.T * wmat_full * N

    sobolev_control = generalized_row_constant(
        residual_rows,
        high_modes,
        wmat_full,
        mp.mpf(args.tol),
    )
    sobolev_absorb, _svals = positive_inverse_constant(
        high_a_modes.T * wmat_n * high_a_modes,
        high_a_modes.T * A * high_a_modes,
        mp.mpf(args.tol),
    )
    source_absorb, source_size = direct_absorption_constant(
        residual_on_n,
        high_a_modes,
        high_avals,
    )
    source_absorb_full, source_size_full = direct_absorption_constant(
        residual_on_n,
        a_modes,
        avals,
    )
    source_absorb_frac = source_absorb / source_absorb_full if source_absorb_full else mp.mpf("0")

    per_source = []
    for t, pstar, rows in blocks:
        row_n = rows * N
        eta, size = direct_absorption_constant(row_n, high_a_modes, high_avals)
        full_eta, _full_size = direct_absorption_constant(row_n, a_modes, avals)
        hconst, _ = graph_constant(rows, high_modes, wmat_full, mp.mpf(args.tol))
        per_source.append(
            {
                "t": f(t),
                "pStarH": f(pstar),
                "sobolevControl": f(hconst),
                "sourceAbsorb": f(eta),
                "sourceAbsorbFull": f(full_eta),
                "sourceAbsorbFrac": f(eta / full_eta if full_eta else mp.mpf("0")),
                "sourceSize": f(size),
            }
        )

    product = sobolev_control * sobolev_absorb
    print(
        f"Aux regularizer certificate basis={args.basis} constraints={args.constraints} "
        f"cutoff={cutoff} positive={len(avals)}"
    )
    print(f"  W{args.sobolev_order} controls stacked residual: {fmt(sobolev_control, 10)}")
    print(f"  W{args.sobolev_order} absorption into A:       {fmt(sobolev_absorb, 10)}")
    print(f"  W route product:                 {fmt(product, 10)}")
    print(f"  source-range absorption into A:  {fmt(source_absorb, 10)}")
    print(f"  source-range high/full fraction: {fmt(source_absorb_frac, 10)}")
    print("  per-source absorption:")
    for row in per_source:
        print(
            f"    t={row['t']} W={fmt(mp.mpf(row['sobolevControl']), 8)} "
            f"Ssrc<=eta A eta={fmt(mp.mpf(row['sourceAbsorb']), 8)} "
            f"frac={fmt(mp.mpf(row['sourceAbsorbFrac']), 8)}"
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "sobolevOrder": args.sobolev_order,
        "lambdaMinA": f(avals[0]),
        "lambdaMaxA": f(avals[-1]),
        "sobolevControl": f(sobolev_control),
        "sobolevAbsorb": f(sobolev_absorb),
        "sobolevRouteProduct": f(product),
        "sourceRangeAbsorb": f(source_absorb),
        "sourceRangeAbsorbFull": f(source_absorb_full),
        "sourceRangeAbsorbFrac": f(source_absorb_frac),
        "sourceRangeSize": f(source_size),
        "sourceRangeSizeFull": f(source_size_full),
        "perSource": per_source,
        "interpretation": (
            "Full Sobolev W10 controls the residual but is not absorbable by "
            "the compact positive block. The finite source-range regularizer "
            "is the viable auxiliary object; its absorption constant is the "
            "direct high-block Schur/Douglas constant."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
