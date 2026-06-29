#!/usr/bin/env python3
r"""Check that global interval traces observe the local-tower bad directions.

The local trace tower at s0 leaves source-heavy directions.  The global trace
operator over the interval appears to remove them.  This script makes that
statement quantitative:

1. build the exact local trace-tower constrained high block;
2. form the A-normalized source Gram for the corrected two-row source E_u;
3. extract the top source directions;
4. measure their sampled global interval trace norm.

If the top source directions have large global trace, then the global closed
trace condition is precisely the missing observability input.
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
    trace_matrix,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def constraints_for(args, basis: int, *, local: bool) -> int:
    if local:
        return max(1, basis - args.local_constraint_offset)
    if args.global_constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(args.global_constraint_ratio * basis)))
    if args.global_constraint_rule == "ratio-round":
        return max(1, min(basis - 1, round(args.global_constraint_ratio * basis)))
    if args.global_constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.global_constraint_offset))
    if args.global_constraint_rule == "fixed":
        return args.global_constraints
    raise ValueError(args.global_constraint_rule)


def child_args(args, basis: int, constraints: int) -> SimpleNamespace:
    out = SimpleNamespace(**vars(args))
    out.basis = basis
    out.constraints = constraints
    return out


def source_nodes(args):
    lo = mp.mpf(args.source_min)
    hi = mp.mpf(args.source_max)
    if args.source_grid == 1:
        return [(lo + hi) / 2]
    step = (hi - lo) / (args.source_grid - 1)
    return [lo + i * step for i in range(args.source_grid)]


def diag_scale_cols(mat, scales):
    out = mp.matrix(mat.rows, len(scales))
    for j, scale in enumerate(scales):
        for i in range(mat.rows):
            out[i, j] = mat[i, j] * scale
    return out


def stack_source_rows(args, polys, e_derivs):
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    nodes = source_nodes(args)
    out = mp.matrix(2 * len(nodes), len(polys))
    for block, u in enumerate(nodes):
        rows, _pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        for i in range(2):
            for j in range(len(polys)):
                out[2 * block + i, j] = rows[i, j]
    return nodes, out


def row_norm(vec):
    return mp.sqrt(mp.fsum(vec[i] * vec[i] for i in range(vec.rows)))


def eig_keep(mat, tol_text):
    vals, vecs = mp.eigsy((mat + mat.T) / 2, eigvals_only=False)
    max_abs = max([abs(v) for v in vals] + [mp.mpf("1")])
    tol = mp.mpf(tol_text) * max_abs
    keep = [i for i, val in enumerate(vals) if val > tol]
    return vals, vecs, keep


def top_indices(vals, count):
    return sorted(range(len(vals)), key=lambda i: vals[i], reverse=True)[:count]


def compute_case(args, basis: int, e_derivs):
    local_constraints = constraints_for(args, basis, local=True)
    global_constraints = constraints_for(args, basis, local=False)
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
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]
    if not high_avals:
        return None

    full_high_modes = N * high_a_modes
    scaled_modes = diag_scale_cols(
        full_high_modes,
        [1 / mp.sqrt(lam) for lam in high_avals],
    )

    source_nodes_used, source_rows = stack_source_rows(args, polys, e_derivs)
    source_on_scaled = source_rows * scaled_modes
    source_gram = source_on_scaled.T * source_on_scaled
    svals, svecs, _skeep = eig_keep(source_gram, args.eig_tol)

    gargs = make_qargs(child_args(args, basis, global_constraints))
    centers, R_global = trace_matrix(polys, gargs)
    global_on_scaled = R_global * scaled_modes
    trace_gram = global_on_scaled.T * global_on_scaled
    tvals, _tvecs, tkeep = eig_keep(trace_gram, args.eig_tol)

    top_rows = []
    for idx in top_indices(svals, min(args.top, len(svals))):
        z = mp.matrix(svecs.rows, 1)
        for i in range(svecs.rows):
            z[i] = svecs[i, idx]
        trace_vec = global_on_scaled * z
        source_vec = source_on_scaled * z
        trace_l2 = row_norm(trace_vec)
        source_l2 = row_norm(source_vec)
        trace_max = max(abs(trace_vec[i]) for i in range(trace_vec.rows)) if trace_vec.rows else mp.mpf("0")
        top_rows.append(
            {
                "sourceMode": idx,
                "sourceEigenvalue": f(svals[idx]),
                "sourceNorm": f(source_l2),
                "globalTraceL2": f(trace_l2),
                "globalTraceMax": f(trace_max),
                "tracePerSourceNorm": f(trace_l2 / source_l2 if source_l2 else mp.mpf("0")),
            }
        )

    kernel_idx = [i for i, val in enumerate(tvals) if val <= mp.mpf(args.kernel_tol)]
    kernel_source_max = mp.mpf("0")
    if kernel_idx:
        kernel_basis = columns(_tvecs, kernel_idx)
        restricted = kernel_basis.T * source_gram * kernel_basis
        kvals = mp.eigsy((restricted + restricted.T) / 2, eigvals_only=True)
        kernel_source_max = max(mp.mpf("0"), kvals[-1])

    return {
        "basis": basis,
        "localTowerConstraints": local_constraints,
        "globalTraceConstraints": global_constraints,
        "globalTraceMin": f(centers[0]) if centers else None,
        "globalTraceMax": f(centers[-1]) if centers else None,
        "localRank": local_rank,
        "localNullity": local_nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinHigh": f(high_avals[0]),
        "lambdaMaxHigh": f(high_avals[-1]),
        "sourceTopEigenvalue": f(max(svals) if len(svals) else mp.mpf("0")),
        "globalTraceRankOnLocalHigh": len(tkeep),
        "globalTraceKernelDimOnLocalHigh": len(kernel_idx),
        "globalTraceMinPositiveEigenvalue": f(min((tvals[i] for i in tkeep), default=mp.mpf("0"))),
        "globalTraceMaxEigenvalue": f(max(tvals) if len(tvals) else mp.mpf("0")),
        "sourceConstantOnGlobalTraceKernel": f(kernel_source_max),
        "topSourceDirections": top_rows,
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
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--top", type=int, default=4)
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
    parser.add_argument("--eig-tol", default="1e-30")
    parser.add_argument("--kernel-tol", default="1e-24")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="global_trace_kills_local_bad_modes.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_constraints = max(constraints_for(args, basis, local=True) for basis in bases)
    max_q = max(args.max_trace_q, max_constraints - 1)
    print(
        f"Global trace kills local bad modes bases={bases} source_grid={args.source_grid}",
        flush=True,
    )
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis local global high src_top trace_rank trace_ker top_trace/source", flush=True)

    rows = []
    for basis in bases:
        row = compute_case(args, basis, e_derivs)
        if row is None:
            continue
        rows.append(row)
        top_ratio = row["topSourceDirections"][0]["tracePerSourceNorm"] if row["topSourceDirections"] else 0.0
        print(
            f"  {basis:5d} {row['localTowerConstraints']:5d} "
            f"{row['globalTraceConstraints']:6d} {row['highModes']:4d} "
            f"{fmt(mp.mpf(row['sourceTopEigenvalue']), 8):>11} "
            f"{row['globalTraceRankOnLocalHigh']:10d} "
            f"{row['globalTraceKernelDimOnLocalHigh']:9d} "
            f"{fmt(mp.mpf(top_ratio), 8):>16}",
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
        "localConstraintOffset": args.local_constraint_offset,
        "globalConstraintRule": args.global_constraint_rule,
        "globalConstraintRatio": args.global_constraint_ratio,
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
            "Top source-heavy directions left by the local trace tower are "
            "measured by the global interval trace operator.  Large trace "
            "norms and zero/small trace-kernel source constant support the "
            "global observability theorem."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
