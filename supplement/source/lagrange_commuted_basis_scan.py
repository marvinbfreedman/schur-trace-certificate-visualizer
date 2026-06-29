#!/usr/bin/env python3
r"""Basis-refinement scan for the commuted Sturm domination certificate.

The previous finite theorem closed at one Galerkin size:

    ||E_u f||^2 <= H_{m,M}(u) ||f||_{W^{m,2}}^2,
    ||f||_{W^{m,2}}^2 <= C_{m,M} <S_m^comm f,f>.

This script recomputes both constants as the polynomial basis changes and
reports

    D_{m,M} = max_u H_{m,M}(u) C_{m,M}.

Stability of D_{m,M}<1 after a fixed cutoff M is the finite evidence needed
before trying to prove the continuum commuted Sturm elliptic estimate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_commuted_kernel_energy import (  # noqa: E402
    commuted_kernel_matrix,
    derivative_matrix,
    positive_inverse_constant,
)
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import (  # noqa: E402
    graph_constant,
    residual_rows_for,
    sobolev_matrix,
)
from lagrange_split_control_certificate import parse_ints  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    trace_matrix,
)
from trace_lagrange_adjoint_control import load_exact_trace  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def child_args(args, basis: int) -> SimpleNamespace:
    out = SimpleNamespace(**vars(args))
    out.basis = basis
    if args.constraints_mode == "ratio":
        out.constraints = max(1, int(mp.floor(mp.mpf(args.constraint_ratio) * basis)))
    else:
        out.constraints = int(args.constraints)
    return out


def max_hardy_constant(args, polys, modes, cutoff: int, order: int):
    graph = sobolev_matrix(polys, mp.mpf(args.L), order, top_only=False)
    high_modes = columns(modes, list(range(cutoff, modes.cols)))
    if high_modes.cols == 0:
        return mp.mpf("0"), None

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    out = mp.mpf("0")
    worst_t = None
    for text in args.t_values.replace(",", " ").split():
        t = mp.mpf(text)
        residual_rows, _pstar = residual_rows_for(
            args,
            polys,
            e_derivs,
            r_nodes,
            r_weights,
            t,
        )
        const, _graph_max = graph_constant(
            residual_rows,
            high_modes,
            graph,
            mp.mpf(args.graph_tol),
        )
        if const > out:
            out = const
            worst_t = t
    return out, worst_t


def commuted_graph_constant(args, K, polys, N, a_modes, cutoff: int, order: int):
    high = columns(a_modes, list(range(cutoff, a_modes.cols)))
    if high.cols == 0:
        return mp.mpf("0"), mp.mpf("0"), mp.mpf("0")
    W = N.T * sobolev_matrix(polys, mp.mpf(args.L), order, top_only=False) * N
    D = derivative_matrix(polys, mp.mpf(args.L))
    Sfull = commuted_kernel_matrix(K, D, order)
    S = N.T * Sfull * N
    Wh = high.T * W * high
    Sh = high.T * S * high
    const, svals = positive_inverse_constant(Wh, Sh, mp.mpf(args.denom_tol))
    return const, svals[0], svals[-1]


def compute_basis(args, basis: int):
    bargs = child_args(args, basis)
    qargs = make_qargs(bargs)
    K, polys = gram_matrix(qargs, mp.mpf(bargs.omega), mp.mpf(bargs.L))
    _centers, R = trace_matrix(polys, qargs)
    N, _U, rank, nullity = split_spaces(R, bargs.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(bargs.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    modes = N * a_modes

    orders = parse_ints(bargs.orders)
    cutoffs = parse_ints(bargs.cutoffs)
    graph_mats = {
        order: sobolev_matrix(polys, mp.mpf(bargs.L), order, top_only=False)
        for order in orders
    }
    quotient_graph_mats = {
        order: N.T * graph_mats[order] * N
        for order in orders
    }
    D = derivative_matrix(polys, mp.mpf(bargs.L))
    commuted_mats = {
        order: N.T * commuted_kernel_matrix(K, D, order) * N
        for order in orders
    }
    _vals, e_derivs, _lam_derivs = load_exact_trace(bargs)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(bargs.kernel_rmax), bargs.kernel_order)
    residual_by_t = []
    for text in bargs.t_values.replace(",", " ").split():
        t = mp.mpf(text)
        residual_rows, _pstar = residual_rows_for(
            bargs,
            polys,
            e_derivs,
            r_nodes,
            r_weights,
            t,
        )
        residual_by_t.append((t, residual_rows))

    rows = []
    for cutoff in cutoffs:
        cutoff = min(cutoff, len(avals))
        high_modes = columns(modes, list(range(cutoff, modes.cols)))
        high_a_modes = columns(a_modes, list(range(cutoff, a_modes.cols)))
        for order in orders:
            hconst = mp.mpf("0")
            worst_t = None
            if high_modes.cols:
                for t, residual_rows in residual_by_t:
                    const, _graph_max = graph_constant(
                        residual_rows,
                        high_modes,
                        graph_mats[order],
                        mp.mpf(bargs.graph_tol),
                    )
                    if const > hconst:
                        hconst = const
                        worst_t = t

            cconst = mp.mpf("0")
            smin = mp.mpf("0")
            smax = mp.mpf("0")
            if high_a_modes.cols:
                Wh = high_a_modes.T * quotient_graph_mats[order] * high_a_modes
                Sh = high_a_modes.T * commuted_mats[order] * high_a_modes
                cconst, svals = positive_inverse_constant(
                    Wh,
                    Sh,
                    mp.mpf(bargs.denom_tol),
                )
                smin = svals[0]
                smax = svals[-1]
            product = hconst * cconst
            rows.append(
                {
                    "cutoff": cutoff,
                    "order": order,
                    "hardyMax": f(hconst),
                    "hardyWorstT": f(worst_t) if worst_t is not None else None,
                    "commutedConstant": f(cconst),
                    "product": f(product),
                    "normProduct": f(mp.sqrt(max(mp.mpf("0"), product))),
                    "sMin": f(smin),
                    "sMax": f(smax),
                    "subunit": bool(product < 1),
                    "emptyHighBlock": cutoff >= len(avals),
                }
            )

    return {
        "basis": basis,
        "constraints": bargs.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "lambdaMinA": f(avals[0]) if avals else None,
        "lambdaMaxA": f(avals[-1]) if avals else None,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="16 18 20")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=11)
    parser.add_argument("--constraints-mode", choices=["ratio", "fixed"], default="ratio")
    parser.add_argument("--constraint-ratio", default="0.625")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoffs", default="4 5")
    parser.add_argument("--orders", default="8 10")
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
    parser.add_argument("--graph-tol", default="1e-40")
    parser.add_argument("--denom-tol", default="1e-40")
    parser.add_argument("--json-out", default="lagrange_commuted_basis_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    results = []

    print("Commuted domination basis scan", flush=True)
    print("  basis cons null pos cutoff order  product   sqrt(product) subunit", flush=True)
    for basis in bases:
        result = compute_basis(args, basis)
        results.append(result)
        for row in result["rows"]:
            print(
                f"  {basis:5d} {result['constraints']:4d} {result['nullity']:4d} "
                f"{result['positiveModes']:3d} {row['cutoff']:6d} {row['order']:5d} "
                f"{fmt(mp.mpf(row['product']), 8):>10} "
                f"{fmt(mp.mpf(row['normProduct']), 8):>13} "
                f"{'yes' if row['subunit'] else 'no'}",
                flush=True,
            )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "bases": bases,
        "constraintMode": args.constraints_mode,
        "constraintRatio": f(mp.mpf(args.constraint_ratio)),
        "cutoffs": parse_ints(args.cutoffs),
        "orders": parse_ints(args.orders),
        "results": results,
        "interpretation": (
            "Basis-refinement scan for D_{m,M}=max_u H_{m,M}(u) C_{m,M}. "
            "Subunit rows are finite high-block domination certificates."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
