#!/usr/bin/env python3
r"""Hardy/Sobolev graph-norm certificate for the Lagrange residual tail.

After the low Schur block is removed, the remaining theorem should be a trace
estimate of Hardy/Sobolev type.  For the two residual rows

    E_u f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0)),

this script computes the sharp finite constants

    ||E_u f||^2 <= C_{m,M}(u) ||f||_{W^{m,2}(0,L)}^2,
        f in H_M,

where H_M is the span of A-eigenmodes after the first M positive modes and

    ||f||_{W^{m,2}}^2 = sum_{r=0}^m int_0^L |f^(r)(s)|^2 ds.

This is not the final continuum Hardy theorem, but it tests the exact finite
mechanism: the residual rows are local jets, so they should be controlled by a
fixed Sobolev/commuted graph norm on the high block.
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
from lagrange_split_control_certificate import parse_ints  # noqa: E402
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


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def poly_derivative_coeffs(poly, deriv: int):
    if deriv == 0:
        return list(poly)
    if deriv >= len(poly):
        return [mp.mpf("0")]
    return [
        poly[p] * mp.factorial(p) / mp.factorial(p - deriv)
        for p in range(deriv, len(poly))
    ]


def poly_inner(a, b, length):
    total = mp.mpf("0")
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            total += ai * bj * (length ** (i + j + 1)) / (i + j + 1)
    return total


def sobolev_matrix(polys, length, order: int, top_only: bool = False):
    out = mp.matrix(len(polys))
    derivs = [order] if top_only else list(range(order + 1))
    for r in derivs:
        dpolys = [poly_derivative_coeffs(poly, r) for poly in polys]
        for i in range(len(polys)):
            for j in range(i + 1):
                val = poly_inner(dpolys[i], dpolys[j], length)
                out[i, j] += val
                if i != j:
                    out[j, i] += val
    return (out + out.T) / 2


def positive_inverse(mat, tol):
    vals, vecs = mp.eigsy((mat + mat.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > tol]
    inv = mp.matrix(mat.rows)
    if keep:
        vkeep = columns(vecs, keep)
        diag = mp.matrix(len(keep))
        for j, idx in enumerate(keep):
            diag[j, j] = 1 / vals[idx]
        inv = vkeep * diag * vkeep.T
    return vals, keep, inv


def residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, t):
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
    out = mp.matrix(2, len(polys))
    for j in range(len(polys)):
        out[0, j] = brow_poly[0, j]
        out[1, j] = eval_poly[0, j]
    return out, pstar


def graph_constant(residual_rows, high_modes, graph_matrix, tol):
    if high_modes.cols == 0:
        return mp.mpf("0"), mp.mpf("0")
    coeffs = residual_rows * high_modes
    graph = high_modes.T * graph_matrix * high_modes
    _vals, _keep, graph_plus = positive_inverse(graph, tol)
    control = coeffs * graph_plus * coeffs.T
    return max_eig_or_zero(control), max_eig_or_zero(graph)


def residual_a_tail(residual_rows, high_modes, avals):
    if high_modes.cols == 0:
        return mp.mpf("0")
    coeffs = residual_rows * high_modes
    out = mp.matrix(2)
    for j, lam in enumerate(avals):
        for a in range(2):
            for b in range(2):
                out[a, b] += coeffs[a, j] * coeffs[b, j] / lam
    return max_eig_or_zero((out + out.T) / 2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=11)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoffs", default="2 3 4 5")
    parser.add_argument("--sobolev-orders", default="2 4 6 8 10")
    parser.add_argument("--top-only", action="store_true")
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
    parser.add_argument("--json-out", default="lagrange_hardy_graph_certificate.json")
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
    modes = N * columns(avecs_all, keep)

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    graph_mats = {
        order: sobolev_matrix(polys, mp.mpf(args.L), order, top_only=args.top_only)
        for order in parse_ints(args.sobolev_orders)
    }
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]

    print(
        f"Lagrange Hardy graph certificate basis={args.basis} constraints={args.constraints} "
        f"nullity={nullity} positive={len(avals)}",
        flush=True,
    )
    print("  cutoff t       A_tail     W2       W4       W6       W8       W10", flush=True)

    rows = []
    for cutoff in parse_ints(args.cutoffs):
        cutoff = min(cutoff, len(avals))
        high_modes = columns(modes, list(range(cutoff, len(avals))))
        high_avals = avals[cutoff:]
        for t in t_values:
            residual_rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, t)
            a_tail = residual_a_tail(residual_rows, high_modes, high_avals)
            graph_rows = []
            for order, gmat in graph_mats.items():
                const, graph_max = graph_constant(
                    residual_rows,
                    high_modes,
                    gmat,
                    mp.mpf(args.graph_tol),
                )
                graph_rows.append(
                    {
                        "order": order,
                        "constant": f(const),
                        "norm": f(mp.sqrt(max(mp.mpf("0"), const))),
                        "graphMaxEigenvalue": f(graph_max),
                    }
                )
            rows.append(
                {
                    "cutoff": cutoff,
                    "t": f(t),
                    "pStarH": f(pstar),
                    "aTailConstant": f(a_tail),
                    "aTailNorm": f(mp.sqrt(max(mp.mpf("0"), a_tail))),
                    "graph": graph_rows,
                }
            )
            by_order = {item["order"]: mp.mpf(item["constant"]) for item in graph_rows}
            print(
                f"  {cutoff:6d} {fmt(t, 5):>6} {fmt(a_tail, 8):>10} "
                f"{fmt(by_order.get(2, mp.mpf('0')), 6):>8} "
                f"{fmt(by_order.get(4, mp.mpf('0')), 6):>8} "
                f"{fmt(by_order.get(6, mp.mpf('0')), 6):>8} "
                f"{fmt(by_order.get(8, mp.mpf('0')), 6):>8} "
                f"{fmt(by_order.get(10, mp.mpf('0')), 6):>8}",
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
        "sobolevOrders": parse_ints(args.sobolev_orders),
        "topOnly": args.top_only,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
