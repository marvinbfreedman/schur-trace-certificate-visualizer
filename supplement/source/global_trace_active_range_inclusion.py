#!/usr/bin/env python3
r"""Construct the finite source-active range inclusion.

The active-gap scan verifies that the sampled interval trace map has no kernel
on the source-active subspace.  In finite dimension this is equivalent to the
range identity

    E_active = C_sample R_global

on that subspace.  This script constructs

    C_sample = E_active (R_active^T R_active)^+ R_active^T

and reports the residual.  The coefficient norms are allowed to be large; the
proof target is exact range inclusion, with quantitative smallness supplied
later by the source-inactive tail estimate.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import (  # noqa: E402
    constraints_for,
    local_case,
    parse_floats,
)
from global_trace_observability_gap import (  # noqa: E402
    child_args,
    f,
    fmt,
    parse_ints,
    stack_source_rows,
    sym,
)
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402


def frob_norm(mat):
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def op_norm(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat.T * mat), eigvals_only=True)
    return mp.sqrt(max(mp.mpf("0"), vals[-1]))


def inverse_spd(mat, tol_text: str):
    vals, vecs = mp.eigsy(sym(mat), eigvals_only=False)
    scale = max([abs(v) for v in vals] + [mp.mpf("1")])
    tol = mp.mpf(tol_text) * scale
    keep = [i for i, val in enumerate(vals) if val > tol]
    inv = mp.matrix(mat.rows)
    if keep:
        vkeep = columns(vecs, keep)
        for col, idx in enumerate(keep):
            for i in range(mat.rows):
                for j in range(mat.rows):
                    inv[i, j] += vkeep[i, col] * vkeep[j, col] / vals[idx]
    return vals, keep, inv, tol


def active_indices(source_vals, source_top, active_tol):
    floor = active_tol * max(mp.mpf("1"), source_top)
    return [i for i, val in enumerate(source_vals) if val > floor], floor


def compute_row(args, base, e_derivs, ratio: float, active_tol: mp.mpf):
    global_constraints = constraints_for(args, base["basis"], ratio, local=False)
    gargs = make_qargs(child_args(args, base["basis"], global_constraints))
    centers, R_global = trace_matrix(base["polys"], gargs)
    active_idx, floor = active_indices(base["sourceVals"], base["sourceTop"], active_tol)
    active_basis = columns(base["sourceVecs"], active_idx)

    source_active = base["sourceOnScaled"] * active_basis
    trace_active = R_global * base["scaledModes"] * active_basis

    gram = trace_active.T * trace_active
    gvals, gkeep, gplus, gtol = inverse_spd(gram, args.trace_tol)
    left_inverse = gplus * trace_active.T
    coeff = source_active * left_inverse
    reconstructed = coeff * trace_active
    residual = source_active - reconstructed
    source_norm = frob_norm(source_active)
    residual_norm = frob_norm(residual)
    coeff_norm = frob_norm(coeff)
    coeff_op = op_norm(coeff)

    return {
        "basis": base["basis"],
        "ratio": ratio,
        "globalTraceConstraints": global_constraints,
        "globalTraceMin": f(centers[0]) if centers else None,
        "globalTraceMax": f(centers[-1]) if centers else None,
        "activeTol": f(active_tol),
        "activeFloor": f(floor),
        "activeDim": len(active_idx),
        "traceRankOnActive": len(gkeep),
        "traceGramMin": f(min(gvals) if len(gvals) else mp.mpf("0")),
        "traceGramMax": f(max(gvals) if len(gvals) else mp.mpf("0")),
        "traceGramTol": f(gtol),
        "sourceActiveFrob": f(source_norm),
        "coefficientFrob": f(coeff_norm),
        "coefficientOp": f(coeff_op),
        "rangeResidualFrob": f(residual_norm),
        "rangeResidualRelative": f(residual_norm / source_norm if source_norm else mp.mpf("0")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18")
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--global-constraint-ratios", default="0.625")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tols", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=70)
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
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="global_trace_active_range_inclusion.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    ratios = parse_floats(args.global_constraint_ratios)
    active_tols = [mp.mpf(piece) for piece in args.active_tols.replace(",", " ").split()]
    max_local_constraints = max(
        constraints_for(args, basis, args.global_constraint_ratio, local=True)
        for basis in bases
    )
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    print(
        f"Global trace active range inclusion bases={bases} ratios={ratios}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis ratio active_tol active_dim rank rel_res coeff_frob", flush=True)

    rows = []
    for basis in bases:
        base = local_case(args, basis, e_derivs)
        _source_nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs)
        base["sourceOnScaled"] = source_rows * base["scaledModes"]
        for ratio in ratios:
            for active_tol in active_tols:
                row = compute_row(args, base, e_derivs, ratio, active_tol)
                rows.append(row)
                print(
                    f"  {basis:5d} {ratio:5.3f} {fmt(active_tol, 4):>10} "
                    f"{row['activeDim']:10d} {row['traceRankOnActive']:4d} "
                    f"{fmt(mp.mpf(row['rangeResidualRelative']), 8):>12} "
                    f"{fmt(mp.mpf(row['coefficientFrob']), 8):>12}",
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
        "globalConstraintRatios": ratios,
        "activeTols": [f(x) for x in active_tols],
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
            "Finite construction of E_active=C_sample R_global on the "
            "source-active block.  The residual should be roundoff when the "
            "trace map has full column rank on the active subspace."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
