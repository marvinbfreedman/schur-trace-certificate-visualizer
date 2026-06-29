#!/usr/bin/env python3
r"""Direct source-row Schur test under the local continuum trace tower.

The individual endpoint jets are too large to control separately in the
compact endpoint-B energy.  The corrected object is the special two-row source
operator

    E_u f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0)).

This script imposes the exact local trace tower

    D_s^q Lambda_s(f)|_{s=s0}=0

and computes the sharp finite constants

    E_u^* E_u <= C(u) A

on the high block after the first M positive A-modes are removed.  It tests
whether the cancellation in the source coefficient family survives the
continuum trace-tower replacement of sampled trace rows.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from aux_regularizer_certificate import direct_absorption_constant  # noqa: E402
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import residual_rows_for  # noqa: E402
from local_trace_tower_representer_scan import (  # noqa: E402
    exact_trace_derivatives,
    tower_matrix,
)
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
)
from source_window_derivative_scan import residual_rows_du_for  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def constraints_for(args, basis: int) -> int:
    if args.constraint_rule == "fixed":
        return args.constraints
    if args.constraint_rule == "basis-minus":
        return max(1, basis - args.constraint_offset)
    if args.constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(args.constraint_ratio * basis)))
    raise ValueError(args.constraint_rule)


def child_args(args, basis: int) -> SimpleNamespace:
    out = SimpleNamespace(**vars(args))
    out.basis = basis
    out.constraints = constraints_for(args, basis)
    return out


def source_nodes(args):
    lo = mp.mpf(args.source_min)
    hi = mp.mpf(args.source_max)
    if args.source_grid == 1:
        return [(lo + hi) / 2]
    step = (hi - lo) / (args.source_grid - 1)
    return [lo + i * step for i in range(args.source_grid)]


def finite_case(args, basis: int, e_derivs):
    bargs = child_args(args, basis)
    qargs = make_qargs(bargs)
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    R = tower_matrix(
        polys,
        mp.mpf(args.s0),
        e_derivs,
        bargs.constraints,
        args.jet_order,
    )
    N, _U, rank, nullity = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    rows = []
    max_high = mp.mpf("0")
    max_d_high = mp.mpf("0")
    max_full = mp.mpf("0")
    max_high_frac = mp.mpf("0")
    max_d_high_frac = mp.mpf("0")
    for u in source_nodes(args):
        e_rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        d_rows, dpstar = residual_rows_du_for(args, polys, e_derivs, r_nodes, r_weights, u)
        e_row_n = e_rows * N
        d_row_n = d_rows * N
        high_eta, high_size = direct_absorption_constant(e_row_n, high_a_modes, high_avals)
        full_eta, full_size = direct_absorption_constant(e_row_n, a_modes, avals)
        d_high_eta, d_high_size = direct_absorption_constant(d_row_n, high_a_modes, high_avals)
        d_full_eta, d_full_size = direct_absorption_constant(d_row_n, a_modes, avals)
        high_frac = high_eta / full_eta if full_eta else mp.mpf("0")
        d_high_frac = d_high_eta / d_full_eta if d_full_eta else mp.mpf("0")
        max_high = max(max_high, high_eta)
        max_d_high = max(max_d_high, d_high_eta)
        max_full = max(max_full, full_eta)
        max_high_frac = max(max_high_frac, high_frac)
        max_d_high_frac = max(max_d_high_frac, d_high_frac)
        rows.append(
            {
                "u": f(u),
                "pStar": f(pstar),
                "dPStar": f(dpstar),
                "highConstant": f(high_eta),
                "fullConstant": f(full_eta),
                "highFullFrac": f(high_frac),
                "dHighConstant": f(d_high_eta),
                "dFullConstant": f(d_full_eta),
                "dHighFullFrac": f(d_high_frac),
                "highRowSize": f(high_size),
                "fullRowSize": f(full_size),
                "dHighRowSize": f(d_high_size),
                "dFullRowSize": f(d_full_size),
            }
        )

    return {
        "basis": basis,
        "towerConstraints": bargs.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinHigh": f(high_avals[0]) if high_avals else None,
        "lambdaMaxHigh": f(high_avals[-1]) if high_avals else None,
        "maxHighConstant": f(max_high),
        "maxDerivativeHighConstant": f(max_d_high),
        "maxFullConstant": f(max_full),
        "maxHighFullFrac": f(max_high_frac),
        "maxDerivativeHighFullFrac": f(max_d_high_frac),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20 22")
    parser.add_argument("--basis", type=int, default=20)
    parser.add_argument("--constraint-rule", choices=["fixed", "basis-minus", "ratio-floor"], default="basis-minus")
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--constraint-offset", type=int, default=8)
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
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
    parser.add_argument("--json-out", default="local_trace_tower_source_row_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_constraints = max(constraints_for(args, basis) for basis in bases)
    max_q = max(args.max_trace_q, max_constraints - 1)

    print(
        f"Local trace tower source-row scan bases={bases} "
        f"source_grid={args.source_grid} max_q={max_q}",
        flush=True,
    )
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis tower pos high max_E_high max_E_frac max_dE_high max_dE_frac", flush=True)

    rows = []
    for basis in bases:
        row = finite_case(args, basis, e_derivs)
        rows.append(row)
        print(
            f"  {basis:5d} {row['towerConstraints']:5d} {row['positiveModes']:3d} "
            f"{row['highModes']:4d} {fmt(mp.mpf(row['maxHighConstant']), 8):>12} "
            f"{fmt(mp.mpf(row['maxHighFullFrac']), 8):>10} "
            f"{fmt(mp.mpf(row['maxDerivativeHighConstant']), 8):>12} "
            f"{fmt(mp.mpf(row['maxDerivativeHighFullFrac']), 8):>10}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceGrid": args.source_grid,
        "bases": bases,
        "constraintRule": args.constraint_rule,
        "constraintOffset": args.constraint_offset,
        "cutoff": args.cutoff,
        "jetOrder": args.jet_order,
        "maxTraceQ": max_q,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Direct test of the source-row Schur/observability theorem under "
            "the exact local trace tower.  This preserves the p,b source-row "
            "combination instead of bounding endpoint jets separately."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
