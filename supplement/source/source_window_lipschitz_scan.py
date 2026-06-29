#!/usr/bin/env python3
r"""Operator-Lipschitz scan for the source-window family E_u.

After the source quadrature refinement, the continuum lift should prove that
the two-row source operator

    E_u f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0))

is uniformly controlled on the closed trace high block.  This script checks
the finite-net mechanism in the A-energy metric:

    E_u^* E_u <= eta(u) A,
    (E_v-E_u)^*(E_v-E_u)/(v-u)^2 <= L(u,v)^2 A.

If the second bound is moderate, then a source grid with mesh h gives the
covering estimate

    E_u^*E_u <= 2 max_j eta(u_j) A + 2 (h/2)^2 L^2 A

for every u between neighboring source nodes.  This is the numerical analogue
of the desired continuum source-window Hardy/Green theorem.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=17)
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
    parser.add_argument("--json-out", default="source_window_lipschitz_scan.json")
    args = parser.parse_args()

    if args.source_grid < 2:
        raise SystemExit("--source-grid must be at least 2")

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

    source_min = mp.mpf(args.source_min)
    source_max = mp.mpf(args.source_max)
    step = (source_max - source_min) / (args.source_grid - 1)
    source_nodes = [source_min + i * step for i in range(args.source_grid)]

    sampled = []
    row_mats = []
    for u in source_nodes:
        rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        row_n = rows * N
        high_eta, high_size = direct_absorption_constant(row_n, high_a_modes, high_avals)
        full_eta, full_size = direct_absorption_constant(row_n, a_modes, avals)
        frac = high_eta / full_eta if full_eta else mp.mpf("0")
        sampled.append(
            {
                "u": f(u),
                "pStarH": f(pstar),
                "highAbsorb": f(high_eta),
                "fullAbsorb": f(full_eta),
                "highFullFrac": f(frac),
                "highSize": f(high_size),
                "fullSize": f(full_size),
            }
        )
        row_mats.append(row_n)

    intervals = []
    for i in range(len(source_nodes) - 1):
        du = source_nodes[i + 1] - source_nodes[i]
        diff = (row_mats[i + 1] - row_mats[i]) / du
        high_lip, high_lip_size = direct_absorption_constant(diff, high_a_modes, high_avals)
        full_lip, full_lip_size = direct_absorption_constant(diff, a_modes, avals)
        intervals.append(
            {
                "left": f(source_nodes[i]),
                "right": f(source_nodes[i + 1]),
                "highLipAbsorb": f(high_lip),
                "fullLipAbsorb": f(full_lip),
                "highLipSize": f(high_lip_size),
                "fullLipSize": f(full_lip_size),
            }
        )

    max_high = max(mp.mpf(row["highAbsorb"]) for row in sampled)
    min_full = min(mp.mpf(row["fullAbsorb"]) for row in sampled)
    max_frac = max(mp.mpf(row["highFullFrac"]) for row in sampled)
    max_high_lip = max(mp.mpf(row["highLipAbsorb"]) for row in intervals)
    radius = step / 2
    cover_high = 2 * max_high + 2 * radius * radius * max_high_lip
    cover_frac = cover_high / min_full if min_full else mp.mpf("0")

    print(
        f"Source-window Lipschitz scan basis={args.basis} constraints={args.constraints} "
        f"cutoff={cutoff} grid={args.source_grid}",
        flush=True,
    )
    print(f"  sampled max high/full = {fmt(max_frac, 10)}", flush=True)
    print(f"  max high Lipschitz A-constant = {fmt(max_high_lip, 10)}", flush=True)
    print(f"  mesh radius = {fmt(radius, 8)}", flush=True)
    print(
        f"  covering high/full bound vs min sampled full = {fmt(cover_frac, 10)}",
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
        "cutoff": cutoff,
        "sourceMin": f(source_min),
        "sourceMax": f(source_max),
        "sourceGrid": args.source_grid,
        "mesh": f(step),
        "meshRadius": f(radius),
        "maxSampleHighAbsorb": f(max_high),
        "minSampleFullAbsorb": f(min_full),
        "maxSampleHighFullFrac": f(max_frac),
        "maxHighLipAbsorb": f(max_high_lip),
        "coverHighAbsorb": f(cover_high),
        "coverHighFullFracVsMinFull": f(cover_frac),
        "sampled": sampled,
        "intervals": intervals,
        "interpretation": (
            "Finite-net source-window lift. The derivative row controls the "
            "operator change E_u-E_v in the A-energy metric; the covering "
            "constant is a conservative pointwise continuum bound over the "
            "source window."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
