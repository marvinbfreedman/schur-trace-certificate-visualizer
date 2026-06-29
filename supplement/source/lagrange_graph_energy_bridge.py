#!/usr/bin/env python3
r"""Bridge Sobolev graph control to endpoint energy on the high block.

The finite Hardy graph certificate proves

    ||E_u f||^2 <= C_{m,M}(u) ||f||_{W^{m,2}}^2,      f in H_M.

The next analytic step is to control the Sobolev graph norm by a commuted
Sturm/Volterra energy on H_M.  As a first stress test, this script checks the
naive strongest finite version using the endpoint energy itself:

    ||f||_{W^{m,2}}^2 <= G_{m,M} <A f,f>,             f in H_M.

If G_{m,M} is large, the endpoint energy A alone is not the right graph
energy; the continuum theorem must use the genuinely commuted Sturm energy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import sobolev_matrix  # noqa: E402
from lagrange_split_control_certificate import parse_ints  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, trace_matrix  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def max_generalized_graph_constant(wmat, modes, avals, cutoff):
    if cutoff >= len(avals):
        return mp.mpf("0")
    high = columns(modes, list(range(cutoff, len(avals))))
    graph = high.T * wmat * high
    best = mp.mpf("0")
    for i, lam in enumerate(avals[cutoff:]):
        # In the A-eigenbasis, <A f,f> is diagonal with entries lambda_i.
        for j in range(graph.cols):
            pass
    scaled = mp.matrix(graph.rows)
    for i in range(graph.rows):
        for j in range(graph.cols):
            scaled[i, j] = graph[i, j] / mp.sqrt(avals[cutoff + i] * avals[cutoff + j])
    vals = mp.eigsy((scaled + scaled.T) / 2, eigvals_only=True)
    return max(mp.mpf("0"), vals[-1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=11)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--cutoffs", default="2 3 4 5")
    parser.add_argument("--sobolev-orders", default="2 4 6 8 10")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="lagrange_graph_energy_bridge.json")
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

    rows = []
    print(
        f"Graph-to-energy bridge basis={args.basis} constraints={args.constraints} "
        f"positive={len(avals)}",
        flush=True,
    )
    print("  cutoff       G2          G4          G6          G8          G10", flush=True)
    wmats = {
        order: sobolev_matrix(polys, mp.mpf(args.L), order, top_only=False)
        for order in parse_ints(args.sobolev_orders)
    }
    for cutoff in parse_ints(args.cutoffs):
        row = {"cutoff": cutoff, "constants": []}
        for order, wmat in wmats.items():
            const = max_generalized_graph_constant(wmat, modes, avals, cutoff)
            row["constants"].append(
                {
                    "order": order,
                    "constant": f(const),
                    "norm": f(mp.sqrt(max(mp.mpf("0"), const))),
                }
            )
        rows.append(row)
        by_order = {item["order"]: mp.mpf(item["constant"]) for item in row["constants"]}
        print(
            f"  {cutoff:6d} {fmt(by_order.get(2, 0), 8):>11} "
            f"{fmt(by_order.get(4, 0), 8):>11} {fmt(by_order.get(6, 0), 8):>11} "
            f"{fmt(by_order.get(8, 0), 8):>11} {fmt(by_order.get(10, 0), 8):>11}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "lambdaMin": f(avals[0]),
        "lambdaMax": f(avals[-1]),
        "rows": rows,
        "interpretation": (
            "These constants test whether endpoint energy A alone controls "
            "Sobolev graph norms on H_M. Large constants mean the analytic "
            "proof needs a commuted Sturm graph energy."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
