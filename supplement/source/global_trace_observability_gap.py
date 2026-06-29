#!/usr/bin/env python3
r"""Finite model for the global interval-observability contradiction.

The previous diagnostic showed that local-tower bad directions are measured by
the sampled global interval trace.  This script turns that statement into the
finite algebraic inequality behind the compactness proof.

On the A-normalized high block left by the local trace tower, let

    S = E^*E,        T = R_global^* R_global.

The compactness contradiction asks for the source-active escaping directions to
be observable by the interval trace.  In finite dimension this is:

    S <= C T

modulo the harmless part where both S and T vanish.  Equivalently, the source
constant for the trace-penalized denominator

    I + beta T / ||T||

must collapse as beta grows, with limiting source mass zero on ker(T).

This is still a numerical model, not the continuum proof.  Its purpose is to
identify the exact global trace inequality that the Lagrange/compactness
argument has to prove analytically.
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


def parse_mp_list(text: str) -> list[mp.mpf]:
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


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


def sym(mat):
    return (mat + mat.T) / 2


def max_eig(mat):
    if mat.rows == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max(mp.mpf("0"), vals[-1])


def frob_norm(mat):
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def positive_eig(mat, rel_tol_text):
    vals, vecs = mp.eigsy(sym(mat), eigvals_only=False)
    scale = max([abs(val) for val in vals] + [mp.mpf("1")])
    tol = mp.mpf(rel_tol_text) * scale
    keep = [i for i, val in enumerate(vals) if val > tol]
    zero = [i for i, val in enumerate(vals) if val <= tol]
    return vals, vecs, keep, zero, tol


def generalized_constant(numer, denom, rel_tol_text):
    vals, vecs, keep, zero, tol = positive_eig(denom, rel_tol_text)
    kernel_numer = mp.mpf("0")
    kernel_cross = mp.mpf("0")
    if zero:
        zvecs = columns(vecs, zero)
        kernel_numer = max_eig(zvecs.T * numer * zvecs)
        if keep:
            vkeep = columns(vecs, keep)
            kernel_cross = frob_norm(zvecs.T * numer * vkeep)
    if not keep:
        return mp.inf, vals, kernel_numer, kernel_cross, tol
    vkeep = columns(vecs, keep)
    denom_vals = [vals[i] for i in keep]
    core = vkeep.T * numer * vkeep
    scaled = mp.matrix(len(keep))
    for i in range(len(keep)):
        for j in range(len(keep)):
            scaled[i, j] = core[i, j] / mp.sqrt(denom_vals[i] * denom_vals[j])
    out_vals = mp.eigsy(sym(scaled), eigvals_only=True)
    return max(mp.mpf("0"), out_vals[-1]), vals, kernel_numer, kernel_cross, tol


def penalty_constant(source_gram, trace_normed, beta):
    denom = mp.eye(source_gram.rows) + beta * trace_normed
    const, _vals, _kn, _kc, _tol = generalized_constant(source_gram, denom, "1e-40")
    return const


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
    avals_all, avecs_all = mp.eigsy(sym(A), eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]
    if not high_avals:
        return None

    # Columns are original polynomial coordinates, normalized so the A-energy
    # on this high block is the identity.
    full_high_modes = N * high_a_modes
    scaled_modes = diag_scale_cols(
        full_high_modes,
        [1 / mp.sqrt(lam) for lam in high_avals],
    )

    _source_nodes_used, source_rows = stack_source_rows(args, polys, e_derivs)
    source_on_scaled = source_rows * scaled_modes
    source_gram = sym(source_on_scaled.T * source_on_scaled)
    source_top = max_eig(source_gram)
    svals, svecs = mp.eigsy(source_gram, eigvals_only=False)
    active_floor = mp.mpf(args.source_active_tol) * max(mp.mpf("1"), source_top)
    active_idx = [i for i, val in enumerate(svals) if val > active_floor]
    source_active_basis = columns(svecs, active_idx)

    gargs = make_qargs(child_args(args, basis, global_constraints))
    centers, R_global = trace_matrix(polys, gargs)
    global_on_scaled = R_global * scaled_modes
    trace_gram_raw = sym(global_on_scaled.T * global_on_scaled)
    trace_vals = mp.eigsy(trace_gram_raw, eigvals_only=True)
    trace_top = max(mp.mpf("0"), trace_vals[-1]) if len(trace_vals) else mp.mpf("0")
    trace_normed = trace_gram_raw / trace_top if trace_top else trace_gram_raw

    trace_const_normed, tvals_normed, kernel_source, kernel_cross, trace_tol = generalized_constant(
        source_gram,
        trace_normed,
        args.trace_tol,
    )
    trace_const_raw = trace_const_normed / trace_top if trace_top else mp.inf
    active_trace_const_normed = mp.mpf("0")
    active_kernel_source = mp.mpf("0")
    active_kernel_cross = mp.mpf("0")
    active_trace_min = mp.mpf("0")
    active_trace_max = mp.mpf("0")
    if active_idx:
        active_source = source_active_basis.T * source_gram * source_active_basis
        active_trace = source_active_basis.T * trace_normed * source_active_basis
        active_trace_const_normed, active_tvals, active_kernel_source, active_kernel_cross, _active_tol = generalized_constant(
            active_source,
            active_trace,
            args.trace_tol,
        )
        active_trace_min = min(active_tvals)
        active_trace_max = max(active_tvals)
    beta_rows = []
    for beta in parse_mp_list(args.betas):
        const = penalty_constant(source_gram, trace_normed, beta)
        beta_rows.append(
            {
                "beta": f(beta),
                "constant": f(const),
                "constantFracOfSourceTop": f(const / source_top if source_top else mp.mpf("0")),
            }
        )

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
        "sourceTop": f(source_top),
        "traceTop": f(trace_top),
        "traceMinRaw": f(trace_vals[0]) if len(trace_vals) else None,
        "traceRankTol": f(trace_tol),
        "traceOnlyConstantNormed": f(trace_const_normed),
        "traceOnlyConstantRaw": f(trace_const_raw),
        "sourceOnTraceKernel": f(kernel_source),
        "sourceOnTraceKernelFrac": f(kernel_source / source_top if source_top else mp.mpf("0")),
        "sourceTraceKernelCrossFrob": f(kernel_cross),
        "sourceActiveTol": f(mp.mpf(args.source_active_tol)),
        "sourceActiveFloor": f(active_floor),
        "sourceActiveDim": len(active_idx),
        "sourceActiveEigenvalues": [f(svals[i]) for i in active_idx],
        "activeTraceOnlyConstantNormed": f(active_trace_const_normed),
        "activeSourceOnTraceKernel": f(active_kernel_source),
        "activeSourceOnTraceKernelFrac": f(active_kernel_source / source_top if source_top else mp.mpf("0")),
        "activeSourceTraceKernelCrossFrob": f(active_kernel_cross),
        "activeTraceMinNormed": f(active_trace_min),
        "activeTraceMaxNormed": f(active_trace_max),
        "penalty": beta_rows,
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
    parser.add_argument("--betas", default="0 1e-6 1e-4 1e-2 1 1e2 1e4 1e6 1e8")
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
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--source-active-tol", default="1e-8")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="global_trace_observability_gap.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_constraints = max(constraints_for(args, basis, local=True) for basis in bases)
    max_q = max(args.max_trace_q, max_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    print(
        f"Global trace observability gap bases={bases} source_grid={args.source_grid}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  basis high active source_top active_trace_const active_kernel_frac beta1e6_frac", flush=True)

    rows = []
    for basis in bases:
        row = compute_case(args, basis, e_derivs)
        if row is None:
            continue
        rows.append(row)
        beta1e6 = next(
            (item for item in row["penalty"] if item["beta"] == 1_000_000.0),
            row["penalty"][-1],
        )
        print(
            f"  {basis:5d} {row['highModes']:4d} {row['sourceActiveDim']:6d} "
            f"{fmt(mp.mpf(row['sourceTop']), 8):>12} "
            f"{fmt(mp.mpf(row['activeTraceOnlyConstantNormed']), 8):>18} "
            f"{fmt(mp.mpf(row['activeSourceOnTraceKernelFrac']), 8):>18} "
            f"{fmt(mp.mpf(beta1e6['constantFracOfSourceTop']), 8):>12}",
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
        "betas": [f(x) for x in parse_mp_list(args.betas)],
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "Finite penalty/observability model for the compact global trace "
            "contradiction.  On the A-normalized local bad high block, the "
            "source Gram S is compared with the sampled interval trace Gram T. "
            "Zero source mass on ker T and decay under I+beta T support the "
            "continuum interval-observability theorem."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
