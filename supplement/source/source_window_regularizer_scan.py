#!/usr/bin/env python3
r"""Source-window scan for the high-block residual regularizer.

The finite source-range certificate with four sample points must become a
uniform source-window theorem.  This script replaces the sample set by a
Gauss-Legendre quadrature window in the source parameter u and checks

    int ||E_u f||^2 du <= eta <A f,f>,       f in H_M cap ker R,

as well as the worst single-source estimate

    ||E_u f||^2 <= eta(u) <A f,f>.

The reported "fraction" is the high-block absorption constant divided by the
corresponding full closed-trace Schur budget.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from aux_regularizer_certificate import direct_absorption_constant  # noqa: E402
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import residual_rows_for  # noqa: E402
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


def source_quadrature(lo, hi, count):
    nodes, weights = mp.gauss_quadrature(count, "legendre")
    mid = (lo + hi) / 2
    half = (hi - lo) / 2
    return [mid + half * node for node in nodes], [half * weight for weight in weights]


def stack_weighted(blocks, polys_len):
    out = mp.matrix(2 * len(blocks), polys_len)
    for block_index, (_u, weight, _pstar, rows) in enumerate(blocks):
        root = mp.sqrt(weight)
        for i in range(2):
            for j in range(polys_len):
                out[2 * block_index + i, j] = root * rows[i, j]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-count", type=int, default=9)
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
    parser.add_argument("--json-out", default="source_window_regularizer_scan.json")
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
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    u_nodes, u_weights = source_quadrature(
        mp.mpf(args.source_min),
        mp.mpf(args.source_max),
        args.source_count,
    )

    blocks = []
    per_source = []
    for u, weight in zip(u_nodes, u_weights):
        rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        row_n = rows * N
        high_eta, high_size = direct_absorption_constant(row_n, high_a_modes, high_avals)
        full_eta, full_size = direct_absorption_constant(row_n, a_modes, avals)
        frac = high_eta / full_eta if full_eta else mp.mpf("0")
        per_source.append(
            {
                "u": f(u),
                "weight": f(weight),
                "pStarH": f(pstar),
                "highAbsorb": f(high_eta),
                "fullAbsorb": f(full_eta),
                "highFullFrac": f(frac),
                "highSize": f(high_size),
                "fullSize": f(full_size),
            }
        )
        blocks.append((u, weight, pstar, rows))

    stacked = stack_weighted(blocks, len(polys)) * N
    high_eta, high_size = direct_absorption_constant(stacked, high_a_modes, high_avals)
    full_eta, full_size = direct_absorption_constant(stacked, a_modes, avals)
    frac = high_eta / full_eta if full_eta else mp.mpf("0")

    worst = max(per_source, key=lambda row: mp.mpf(row["highFullFrac"]))
    print(
        f"Source-window regularizer scan basis={args.basis} constraints={args.constraints} "
        f"cutoff={cutoff} sources={args.source_count}"
    )
    print(
        f"  window=[{args.source_min},{args.source_max}] "
        f"integrated high/full={fmt(frac, 10)} high_eta={fmt(high_eta, 10)}"
    )
    print(
        f"  worst single u={fmt(mp.mpf(worst['u']), 8)} "
        f"frac={fmt(mp.mpf(worst['highFullFrac']), 10)} "
        f"high_eta={fmt(mp.mpf(worst['highAbsorb']), 10)}"
    )
    print("  node     high_eta     high/full")
    for row in per_source:
        print(
            f"  {fmt(mp.mpf(row['u']), 7):>7} "
            f"{fmt(mp.mpf(row['highAbsorb']), 8):>12} "
            f"{fmt(mp.mpf(row['highFullFrac']), 8):>11}"
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
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceCount": args.source_count,
        "lambdaMinA": f(avals[0]),
        "lambdaMaxA": f(avals[-1]),
        "integratedHighAbsorb": f(high_eta),
        "integratedFullAbsorb": f(full_eta),
        "integratedHighFullFrac": f(frac),
        "integratedHighSize": f(high_size),
        "integratedFullSize": f(full_size),
        "worstSingle": worst,
        "perSource": per_source,
        "interpretation": (
            "Gauss-window source regularizer. Stability of the high/full "
            "fraction supports replacing sampled source points by a uniform "
            "source-window Hardy/Green estimate."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
