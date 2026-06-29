#!/usr/bin/env python3
r"""Stability scan for the source-active global trace gap.

The expensive quotient in ``global_trace_observability_gap.py`` is not the
right object to tune: its constants are dominated by nearly trace-null,
source-inactive directions.  The qualitative theorem only needs:

    ker(R_global) cap source-active(S) = {0}.

This script scans that statement over active source cutoffs and sampled trace
densities.  It reports the trace Gram spectrum on the source-active subspace,
plus the source mass carried by the full trace kernel.  It deliberately avoids
the ill-conditioned generalized quotient.
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

from global_trace_observability_gap import (  # noqa: E402
    child_args,
    diag_scale_cols,
    f,
    fmt,
    max_eig,
    parse_ints,
    positive_eig,
    stack_source_rows,
    sym,
)
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives, tower_matrix  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, trace_matrix  # noqa: E402


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def constraints_for(args, basis: int, ratio: float, *, local: bool) -> int:
    if local:
        return max(1, basis - args.local_constraint_offset)
    if args.global_constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(ratio * basis)))
    if args.global_constraint_rule == "ratio-round":
        return max(1, min(basis - 1, round(ratio * basis)))
    if args.global_constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.global_constraint_offset))
    if args.global_constraint_rule == "fixed":
        return args.global_constraints
    raise ValueError(args.global_constraint_rule)


def local_case(args, basis: int, e_derivs):
    local_constraints = constraints_for(args, basis, args.global_constraint_ratio, local=True)
    qargs = make_qargs(child_args(args, basis, local_constraints))
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    R_local = tower_matrix(
        polys,
        mp.mpf(args.s0),
        e_derivs,
        local_constraints,
        args.jet_order,
    )
    N, _U, local_rank, local_nullity = split_spaces(R_local, args.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy(sym(A), eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]
    full_high_modes = N * high_a_modes
    scaled_modes = diag_scale_cols(
        full_high_modes,
        [1 / mp.sqrt(lam) for lam in high_avals],
    )
    _source_nodes, source_rows = stack_source_rows(args, polys, e_derivs)
    source_on_scaled = source_rows * scaled_modes
    source_gram = sym(source_on_scaled.T * source_on_scaled)
    source_top = max_eig(source_gram)
    svals, svecs = mp.eigsy(source_gram, eigvals_only=False)
    return {
        "basis": basis,
        "polys": polys,
        "scaledModes": scaled_modes,
        "sourceGram": source_gram,
        "sourceTop": source_top,
        "sourceVals": svals,
        "sourceVecs": svecs,
        "localConstraints": local_constraints,
        "localRank": local_rank,
        "localNullity": local_nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinHigh": high_avals[0] if high_avals else mp.mpf("0"),
        "lambdaMaxHigh": high_avals[-1] if high_avals else mp.mpf("0"),
    }


def source_mass_on_trace_kernel(source_gram, trace_gram, trace_tol_text):
    _tvals, tvecs, _keep, zero, _tol = positive_eig(trace_gram, trace_tol_text)
    if not zero:
        return mp.mpf("0"), 0
    zvecs = columns(tvecs, zero)
    return max_eig(zvecs.T * source_gram * zvecs), len(zero)


def scan_row(args, base, ratio: float, active_tol: mp.mpf):
    global_constraints = constraints_for(args, base["basis"], ratio, local=False)
    gargs = make_qargs(child_args(args, base["basis"], global_constraints))
    centers, R_global = trace_matrix(base["polys"], gargs)
    trace_gram_raw = sym((R_global * base["scaledModes"]).T * (R_global * base["scaledModes"]))
    trace_top = max_eig(trace_gram_raw)
    trace_normed = trace_gram_raw / trace_top if trace_top else trace_gram_raw
    full_kernel_source, full_kernel_dim = source_mass_on_trace_kernel(
        base["sourceGram"],
        trace_normed,
        args.trace_tol,
    )

    active_floor = active_tol * max(mp.mpf("1"), base["sourceTop"])
    active_idx = [i for i, val in enumerate(base["sourceVals"]) if val > active_floor]
    active_trace_min = mp.mpf("0")
    active_trace_max = mp.mpf("0")
    active_kernel_source = mp.mpf("0")
    active_kernel_dim = 0
    if active_idx:
        active_basis = columns(base["sourceVecs"], active_idx)
        active_source = active_basis.T * base["sourceGram"] * active_basis
        active_trace = sym(active_basis.T * trace_normed * active_basis)
        active_tvals = mp.eigsy(active_trace, eigvals_only=True)
        active_trace_min = min(active_tvals)
        active_trace_max = max(active_tvals)
        active_kernel_source, active_kernel_dim = source_mass_on_trace_kernel(
            active_source,
            active_trace,
            args.trace_tol,
        )

    return {
        "basis": base["basis"],
        "ratio": ratio,
        "globalTraceConstraints": global_constraints,
        "globalTraceMin": f(centers[0]) if centers else None,
        "globalTraceMax": f(centers[-1]) if centers else None,
        "activeTol": f(active_tol),
        "activeDim": len(active_idx),
        "activeFloor": f(active_floor),
        "sourceTop": f(base["sourceTop"]),
        "traceTop": f(trace_top),
        "fullTraceKernelDim": full_kernel_dim,
        "fullTraceKernelSource": f(full_kernel_source),
        "fullTraceKernelSourceFrac": f(full_kernel_source / base["sourceTop"] if base["sourceTop"] else mp.mpf("0")),
        "activeTraceMinNormed": f(active_trace_min),
        "activeTraceMaxNormed": f(active_trace_max),
        "activeTraceKernelDim": active_kernel_dim,
        "activeTraceKernelSource": f(active_kernel_source),
        "activeTraceKernelSourceFrac": f(active_kernel_source / base["sourceTop"] if base["sourceTop"] else mp.mpf("0")),
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
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--global-constraint-ratios", default="0.50 0.625 0.75")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tols", default="1e-6 1e-8 1e-10 1e-12")
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
    parser.add_argument("--json-out", default="global_trace_active_gap_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    ratios = parse_floats(args.global_constraint_ratios)
    active_tols = [mp.mpf(piece) for piece in args.active_tols.replace(",", " ").split()]
    max_local_constraints = max(constraints_for(args, basis, args.global_constraint_ratio, local=True) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    print(
        f"Global trace active-gap scan bases={bases} ratios={ratios} active_tols={args.active_tols}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis ratio active_tol active_dim active_min active_ker_frac full_ker_frac", flush=True)

    rows = []
    base_rows = []
    for basis in bases:
        base = local_case(args, basis, e_derivs)
        base_rows.append(
            {
                key: f(value) if isinstance(value, mp.mpf) else value
                for key, value in base.items()
                if key not in {"polys", "scaledModes", "sourceGram", "sourceVals", "sourceVecs"}
            }
        )
        for ratio in ratios:
            for active_tol in active_tols:
                row = scan_row(args, base, ratio, active_tol)
                rows.append(row)
                print(
                    f"  {basis:5d} {ratio:5.3f} {fmt(active_tol, 4):>10} "
                    f"{row['activeDim']:10d} {fmt(mp.mpf(row['activeTraceMinNormed']), 8):>10} "
                    f"{fmt(mp.mpf(row['activeTraceKernelSourceFrac']), 8):>15} "
                    f"{fmt(mp.mpf(row['fullTraceKernelSourceFrac']), 8):>13}",
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
        "baseRows": base_rows,
        "rows": rows,
        "interpretation": (
            "Stability scan for the qualitative interval source-observability "
            "gap.  The target is zero source mass on the trace kernel inside "
            "the source-active subspace, across active cutoffs and trace sample "
            "densities."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
