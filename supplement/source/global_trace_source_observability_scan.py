#!/usr/bin/env python3
r"""Global sampled-trace source-row observability scan.

The local trace tower at s0 is too weak for an absolute bound, even though the
special source rows have tiny high/full Schur fraction.  The corrected theorem
uses the interval trace operator

    R f = (Lambda_a f)_{a in I},        I=[0.02,0.545],

with Galerkin trace samples refining as the basis grows.

This script keeps the corrected two-row source operator

    E_u f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0))

and measures the finite constants

    E_u^* E_u <= C(u) A

on H_M cap ker R_sampled.  It is the direct numerical certificate for the
global interval observability/range theorem.
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
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    trace_matrix,
)
from source_window_derivative_scan import residual_rows_du_for  # noqa: E402
from trace_lagrange_adjoint_control import load_exact_trace  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def constraints_for(args, basis: int) -> int:
    if args.constraint_rule == "fixed":
        return args.constraints
    if args.constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(args.constraint_ratio * basis)))
    if args.constraint_rule == "ratio-round":
        return max(1, min(basis - 1, round(args.constraint_ratio * basis)))
    if args.constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.constraint_offset))
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
    centers, R = trace_matrix(polys, qargs)
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
    max_d_full = mp.mpf("0")
    max_high_frac = mp.mpf("0")
    max_d_high_frac = mp.mpf("0")
    worst_u = None
    worst_du = None
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
        if high_eta > max_high:
            max_high = high_eta
            worst_u = u
        if d_high_eta > max_d_high:
            max_d_high = d_high_eta
            worst_du = u
        max_full = max(max_full, full_eta)
        max_d_full = max(max_d_full, d_full_eta)
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
        "constraints": bargs.constraints,
        "traceMin": f(centers[0]) if centers else None,
        "traceMax": f(centers[-1]) if centers else None,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinA": f(avals[0]) if avals else None,
        "lambdaMaxA": f(avals[-1]) if avals else None,
        "lambdaMinHigh": f(high_avals[0]) if high_avals else None,
        "lambdaMaxHigh": f(high_avals[-1]) if high_avals else None,
        "maxHighConstant": f(max_high),
        "maxDerivativeHighConstant": f(max_d_high),
        "maxFullConstant": f(max_full),
        "maxDerivativeFullConstant": f(max_d_full),
        "maxHighFullFrac": f(max_high_frac),
        "maxDerivativeHighFullFrac": f(max_d_high_frac),
        "worstU": f(worst_u) if worst_u is not None else None,
        "worstDerivativeU": f(worst_du) if worst_du is not None else None,
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
    parser.add_argument("--constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--constraint-offset", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
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
    parser.add_argument("--json-out", default="global_trace_source_observability_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    print(
        f"Global trace source observability scan bases={bases} "
        f"source_grid={args.source_grid} rule={args.constraint_rule}",
        flush=True,
    )
    vals, e_derivs, lam_derivs = load_exact_trace(args)
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis cons pos high max_E_high max_E_frac max_dE_high max_dE_frac", flush=True)

    rows = []
    for basis in bases:
        row = finite_case(args, basis, e_derivs)
        rows.append(row)
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['positiveModes']:3d} "
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
        "constraintRatio": args.constraint_ratio,
        "cutoff": args.cutoff,
        "jetOrder": args.jet_order,
        "maxTraceQ": args.max_trace_q,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Direct source-row Schur/observability scan using the global "
            "sampled interval trace operator.  This is the finite model for "
            "E_u^*E_u <= C(u)A on H_M cap ker R_global."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
