#!/usr/bin/env python3
r"""Fixed jet representers with the local continuum trace tower.

The sampled trace constraints

    Lambda_a(f)=0 at a finite mesh of a-values

are only a Galerkin proxy for the closed continuum trace condition.  Near a
fixed interior point s0, the continuum condition implies the full differential
tower

    D_s^q Lambda_s(f)|_{s=s0}=0,        q=0,1,2,...

where

    Lambda_s(f)=sum_{k=0}^8 e_k(s) f^(k)(s)/k!.

This script replaces sampled trace rows by the exact local tower rows computed
from the confluent eigen-equation for e(s).  It then recomputes the A-Green
representer norms for the fixed endpoint jets f^(j)(s0), j=0..7.  The purpose
is to separate a genuine continuum closed-trace effect from finite sampling
artifacts.
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

from adjoint_eval_representer_certificate import energy_representer  # noqa: E402
from boundary_row_representer_certificate import jet_eval_row, matrix_to_lists  # noqa: E402
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, poly_derivative_value  # noqa: E402
from trace_concomitant_exact_derivatives import (  # noqa: E402
    center_taylor_mats,
    eigen_taylor,
)
from lambda_differential_closure import trace_derivative_row  # noqa: E402
from endpoint_kb_confluent_mp import integrate  # noqa: E402


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


def exact_trace_derivatives(args, max_q: int):
    big_order = args.jet_order + max_q
    big, _segments = integrate(
        "kb",
        mp.pi,
        mp.mpf(args.s0),
        big_order,
        mp.mpf(args.matrix_rmax),
        args.matrix_order,
    )
    mats = center_taylor_mats(big, args.jet_order, max_q)
    vals, e_derivs, lam_derivs = eigen_taylor(mats, max_q)
    return vals, e_derivs, lam_derivs


def tower_row_on_polys(polys, s0, e_derivs, q: int, jet_order: int):
    jet_row = trace_derivative_row(e_derivs, q, jet_order)
    row = mp.matrix(1, len(polys))
    for i, poly in enumerate(polys):
        row[0, i] = mp.fsum(
            jet_row[k] * poly_derivative_value(poly, s0, k)
            for k in range(len(jet_row))
        )
    return row


def tower_matrix(polys, s0, e_derivs, constraints: int, jet_order: int):
    rows = []
    for q in range(constraints):
        rows.append(tower_row_on_polys(polys, s0, e_derivs, q, jet_order))
    out = mp.matrix(constraints, len(polys))
    for i, row in enumerate(rows):
        for j in range(len(polys)):
            out[i, j] = row[0, j]
    return out


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

    jet_records = []
    reps = []
    max_range_defect = mp.mpf("0")
    for deriv in range(args.jet_count):
        row_n = jet_eval_row(polys, mp.mpf(args.s0), deriv) * N
        rep_n, norm2, rel = energy_representer(row_n, A, high_a_modes, high_avals)
        reps.append(rep_n)
        max_range_defect = max(max_range_defect, rel)
        jet_records.append(
            {
                "deriv": deriv,
                "norm2": f(norm2),
                "norm": f(mp.sqrt(max(mp.mpf("0"), norm2))),
                "rangeRelativeDefect": f(rel),
            }
        )

    G = mp.matrix(args.jet_count)
    for i in range(args.jet_count):
        for j in range(args.jet_count):
            G[i, j] = (reps[i].T * A * reps[j])[0, 0]

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
        "maxRangeRelativeDefect": f(max_range_defect),
        "jetRepresenters": jet_records,
        "jetGram": matrix_to_lists(G),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
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
    parser.add_argument("--jet-count", type=int, default=8)
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
    parser.add_argument("--json-out", default="local_trace_tower_representer_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_constraints = max(constraints_for(args, basis) for basis in bases)
    max_q = max(args.max_trace_q, max_constraints - 1)

    print(
        f"Local trace tower representer scan bases={bases} "
        f"s0={args.s0} max_q={max_q}",
        flush=True,
    )
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis tower rank null pos high eval_norm2 j7_norm2 range_def", flush=True)

    rows = []
    for basis in bases:
        row = finite_case(args, basis, e_derivs)
        rows.append(row)
        eval_norm2 = mp.mpf(row["jetRepresenters"][0]["norm2"])
        j7_norm2 = mp.mpf(row["jetRepresenters"][7]["norm2"])
        print(
            f"  {basis:5d} {row['towerConstraints']:5d} {row['rank']:4d} "
            f"{row['nullity']:4d} {row['positiveModes']:3d} {row['highModes']:4d} "
            f"{fmt(eval_norm2, 8):>11} {fmt(j7_norm2, 8):>11} "
            f"{fmt(mp.mpf(row['maxRangeRelativeDefect']), 8):>10}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "bases": bases,
        "constraintRule": args.constraint_rule,
        "constraintOffset": args.constraint_offset,
        "cutoff": args.cutoff,
        "jetOrder": args.jet_order,
        "jetCount": args.jet_count,
        "maxTraceQ": max_q,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "The sampled trace mesh has been replaced by the exact local "
            "continuum trace tower D_s^q Lambda_s0.  Stabilization here is "
            "evidence for the closed-trace endpoint RKHS/Sobolev representer "
            "bound; growth would indicate that additional Schur/observability "
            "input is still missing."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
