#!/usr/bin/env python3
r"""Basis-refinement scan for the Lagrange low/high split.

This script reruns the residual split certificate while increasing the
Galerkin basis and the sampled trace constraints together.  The goal is to
check whether the low Schur block size is stable: after removing a small fixed
number of low A-modes, the residual high tail should remain tiny as the basis
grows.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_split_control_certificate import (  # noqa: E402
    compute_case,
    parse_floats,
    parse_ints,
)
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


def constraints_for(args, basis: int):
    if args.constraint_rule == "ratio":
        return max(1, min(basis - 1, int(round(args.constraint_ratio * basis))))
    if args.constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.constraint_offset))
    if args.constraint_rule == "target-nullity":
        return max(1, min(basis - 1, basis - args.target_nullity))
    raise ValueError(args.constraint_rule)


def case_args(base, basis: int):
    constraints = constraints_for(base, basis)
    return SimpleNamespace(
        s0=base.s0,
        t_values=base.t_values,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=base.quad_base + base.quad_step * basis,
        laguerre=base.laguerre,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        max_trace_q=base.max_trace_q,
        cutoffs=base.cutoffs,
        moments=base.moments,
        dps=base.dps,
        matrix_order=base.matrix_order,
        matrix_rmax=base.matrix_rmax,
        kernel_order=base.kernel_order,
        kernel_rmax=base.kernel_rmax,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
    )


def summarize_cutoff(row, cutoff: int):
    vals = []
    for source in row["rows"]:
        match = next((item for item in source["splits"] if item["cutoff"] == cutoff), None)
        if match is not None:
            vals.append(match["highFrac"])
    if not vals:
        return None
    return {
        "cutoff": cutoff,
        "minHighFrac": f(min(mp.mpf(v) for v in vals)),
        "maxHighFrac": f(max(mp.mpf(v) for v in vals)),
    }


def first_cutoff_below(row, threshold):
    threshold = mp.mpf(threshold)
    for cutoff in parse_ints(row["cutoffsText"]):
        summary = summarize_cutoff(row, cutoff)
        if summary is None:
            continue
        if mp.mpf(summary["maxHighFrac"]) <= threshold:
            return cutoff
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="12 14 16 18")
    parser.add_argument("--constraint-rule", choices=["ratio", "offset", "target-nullity"], default="ratio")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--constraint-offset", type=int, default=6)
    parser.add_argument("--target-nullity", type=int, default=6)
    parser.add_argument("--quad-base", type=int, default=6)
    parser.add_argument("--quad-step", type=int, default=1)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoffs", default="1 2 3 4 5 6")
    parser.add_argument("--moments", default="0 1 2")
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
    parser.add_argument("--thresholds", default="1e-2 1e-3 1e-4")
    parser.add_argument("--json-out", default="lagrange_split_basis_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]

    rows = []
    print(
        f"Lagrange split basis scan rule={args.constraint_rule} "
        f"basis={args.basis}",
        flush=True,
    )
    print("  basis cons null pos lambda_min  cut<=1e-2 cut<=1e-3 cut<=1e-4 high@4 high@5", flush=True)

    for basis in parse_ints(args.basis):
        cargs = case_args(args, basis)
        qargs = make_qargs(cargs)
        K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
        _centers, R = trace_matrix(polys, qargs)
        N, _U, rank, nullity = split_spaces(R, args.rank_tol)
        A = N.T * K * N
        avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
        keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
        if not keep:
            continue
        avals = [avals_all[i] for i in keep]
        modes = columns(avecs_all, keep)

        source_rows = [
            compute_case(cargs, t, polys, N, modes, avals, e_derivs, r_nodes, r_weights)
            for t in t_values
        ]
        row = {
            "basis": basis,
            "constraints": cargs.constraints,
            "rank": rank,
            "nullity": nullity,
            "positiveModes": len(avals),
            "lambdaMin": f(avals[0]),
            "lambdaMax": f(avals[-1]),
            "cutoffsText": args.cutoffs,
            "cutoffSummaries": [
                summarize_cutoff({"rows": source_rows, "cutoffsText": args.cutoffs}, cutoff)
                for cutoff in parse_ints(args.cutoffs)
            ],
            "thresholdCutoffs": {
                threshold: first_cutoff_below(
                    {"rows": source_rows, "cutoffsText": args.cutoffs},
                    threshold,
                )
                for threshold in args.thresholds.replace(",", " ").split()
            },
            "rows": source_rows,
        }
        rows.append(row)
        summaries = {
            item["cutoff"]: item for item in row["cutoffSummaries"] if item is not None
        }
        high4 = summaries.get(4, {}).get("maxHighFrac")
        high5 = summaries.get(5, {}).get("maxHighFrac")
        print(
            f"  {basis:5d} {cargs.constraints:4d} {nullity:4d} {len(avals):3d} "
            f"{fmt(avals[0], 6):>10} "
            f"{str(row['thresholdCutoffs'].get('1e-2')):>9} "
            f"{str(row['thresholdCutoffs'].get('1e-3')):>9} "
            f"{str(row['thresholdCutoffs'].get('1e-4')):>9} "
            f"{fmt(mp.mpf(high4), 6) if high4 else 'n/a':>8} "
            f"{fmt(mp.mpf(high5), 6) if high5 else 'n/a':>8}",
            flush=True,
        )

    data = {
        "s0": f(mp.mpf(args.s0)),
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "cutoffs": parse_ints(args.cutoffs),
        "moments": parse_floats(args.moments),
        "thresholds": parse_floats(args.thresholds),
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
