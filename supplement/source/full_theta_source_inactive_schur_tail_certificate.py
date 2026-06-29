#!/usr/bin/env python3
r"""Full-theta source-inactive Schur-tail certificate.

This is the finite/quantitative model for the remaining theorem

    ||(I-P_delta) E_Phi f||^2 <= epsilon_delta <A f,f>,
    f in H_M cap ker R_global.

The script works on the same A-normalized high block used by the source-side
rank theorem.  It builds the literal finite full-theta source Gram for
Psi_8=Phi_{<=8}, splits off the top two source-active eigenmodes, and measures
the source-inactive part against the sampled global trace kernel.

The continuum/full-Phi promotion uses the certified source quadrature/tail
operator error from ``full_theta_source_quadrature_consequence_theorem.json``:

    lambda_3(S_Phi) <= lambda_3(S_{8,h}) + eps.

Thus the source-inactive continuum bound is an honest perturbative upper bound,
while the comparison with the finite Schur-tail budget remains diagnostic until
the analytic high-frequency Hardy/Schur theorem is proved.
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

from determinant_gap_bound_diagnostic import sym  # noqa: E402
from full_theta_source_tail_certificate import (  # noqa: E402
    full_weights,
    pieces_from_mode_weights,
    weighted_source_matrix_for_pieces,
)
from global_trace_active_gap_scan import constraints_for, local_case, source_mass_on_trace_kernel  # noqa: E402
from global_trace_observability_gap import child_args, f, fmt  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, diag, trace_matrix  # noqa: E402


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def op_norm_sym(mat):
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def eigensplit_source(source_gram, active_count: int):
    vals, vecs = mp.eigsy(sym(source_gram), eigvals_only=False)
    if len(vals) <= active_count:
        raise ValueError("not enough source eigenvalues for active/inactive split")
    inactive_idx = list(range(len(vals) - active_count))
    active_idx = list(range(len(vals) - active_count, len(vals)))
    inactive_basis = columns(vecs, inactive_idx)
    active_basis = columns(vecs, active_idx)
    inactive_diag = diag([vals[i] for i in inactive_idx])
    active_diag = diag([vals[i] for i in active_idx])
    inactive_gram = inactive_basis * inactive_diag * inactive_basis.T
    active_gram = active_basis * active_diag * active_basis.T
    return vals, vecs, inactive_idx, active_idx, inactive_gram, active_gram


def bridge_tail_budget(path: str | None):
    if not path:
        return None
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in (
        "finiteSchurTailBudgetDiagnostic",
        "sourceInactiveTailBudgetDiagnostic",
    ):
        value = data.get(key)
        if value:
            return value
    for item in data.get("remainingAnalyticItems", []):
        if item.get("name") == "source-inactive high-frequency Schur tail":
            tail = item.get("finiteDiagnostic", {}).get("bestObservedNonzeroOperatorTail")
            if tail:
                return tail
    return None


def compute_case(args, basis: int, ratio: float, e_derivs):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    bargs.source_grid = min(args.source_grid, args.base_source_grid)
    base = local_case(bargs, basis, e_derivs)

    pieces8 = pieces_from_mode_weights(full_weights(args.full_nmax))
    nodes, weights, source_on_scaled, source_gram = weighted_source_matrix_for_pieces(
        args,
        base,
        e_derivs,
        pieces8,
    )
    source_norm = op_norm_sym(source_gram)
    vals, vecs, inactive_idx, active_idx, inactive_gram, active_gram = eigensplit_source(
        source_gram,
        args.active_count,
    )
    active_min = vals[active_idx[0]]
    complement_top = vals[inactive_idx[-1]]
    spectral_gap = active_min - complement_top

    global_constraints = constraints_for(args, basis, ratio, local=False)
    gargs = make_qargs(child_args(args, basis, global_constraints))
    centers, R_global = trace_matrix(base["polys"], gargs)
    trace_on_scaled = R_global * base["scaledModes"]
    trace_gram_raw = sym(trace_on_scaled.T * trace_on_scaled)
    trace_top = op_norm_sym(trace_gram_raw)
    trace_normed = trace_gram_raw / trace_top if trace_top else trace_gram_raw

    active_basis = columns(vecs, active_idx)
    active_trace = sym(active_basis.T * trace_normed * active_basis)
    active_trace_vals = mp.eigsy(active_trace, eigvals_only=True)
    active_source = sym(active_basis.T * source_gram * active_basis)
    active_kernel_source, active_kernel_dim = source_mass_on_trace_kernel(
        active_source,
        active_trace,
        args.trace_tol,
    )
    full_kernel_source, full_kernel_dim = source_mass_on_trace_kernel(
        source_gram,
        trace_normed,
        args.trace_tol,
    )
    inactive_kernel_source, inactive_kernel_dim = source_mass_on_trace_kernel(
        inactive_gram,
        trace_normed,
        args.trace_tol,
    )

    eps = mp.mpf(args.continuum_error)
    continuum_top_lower = max(mp.mpf("0"), source_norm - eps)
    continuum_inactive_upper = complement_top + eps
    continuum_inactive_frac_upper = (
        continuum_inactive_upper / continuum_top_lower if continuum_top_lower else mp.inf
    )

    return {
        "basis": basis,
        "ratio": ratio,
        "globalTraceConstraints": global_constraints,
        "globalTraceMin": f(centers[0]) if centers else None,
        "globalTraceMax": f(centers[-1]) if centers else None,
        "highModes": base["scaledModes"].cols,
        "sourceGrid": args.source_grid,
        "sourceNodeFirstLast": [f(nodes[0]), f(nodes[-1])],
        "sourceTop": f(source_norm),
        "activeEigenvalues": [f(vals[i]) for i in active_idx],
        "complementTopEigenvalue": f(complement_top),
        "spectralGapToActive": f(spectral_gap),
        "finiteInactiveFracOfTop": f(complement_top / source_norm if source_norm else mp.inf),
        "continuumErrorBound": f(eps),
        "continuumTopLower": f(continuum_top_lower),
        "continuumInactiveUpper": f(continuum_inactive_upper),
        "continuumInactiveFracUpper": f(continuum_inactive_frac_upper),
        "traceTop": f(trace_top),
        "activeTraceMinNormed": f(min(active_trace_vals) if len(active_trace_vals) else mp.mpf("0")),
        "activeTraceKernelDim": active_kernel_dim,
        "activeTraceKernelSourceFrac": f(active_kernel_source / source_norm if source_norm else mp.inf),
        "fullTraceKernelDim": full_kernel_dim,
        "fullTraceKernelSourceFrac": f(full_kernel_source / source_norm if source_norm else mp.inf),
        "inactiveTraceKernelDim": inactive_kernel_dim,
        "inactiveTraceKernelSourceFrac": f(inactive_kernel_source / source_norm if source_norm else mp.inf),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--bases", default="22")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=257)
    parser.add_argument("--base-source-grid", type=int, default=9)
    parser.add_argument("--full-nmax", type=int, default=8)
    parser.add_argument("--source-order", type=int, default=10)
    parser.add_argument("--source-rmax", default="8")
    parser.add_argument("--active-count", type=int, default=2)
    parser.add_argument("--global-constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--global-constraint-ratios", default="0.625 0.75")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--dps", type=int, default=55)
    parser.add_argument("--matrix-order", type=int, default=40)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=16)
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
    parser.add_argument("--quadrature-json", default="full_theta_source_quadrature_consequence_theorem.json")
    parser.add_argument("--bridge-json", default="global_schur_tail_budget_consequence_theorem.json")
    parser.add_argument("--json-out", default="full_theta_source_inactive_schur_tail_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    qdata = json.loads(Path(args.quadrature_json).read_text(encoding="utf-8"))
    args.continuum_error = str(qdata["totalContinuumErrorBound"])
    if int(qdata["sourceGrid"]) != args.source_grid:
        print(
            f"warning: quadrature JSON grid {qdata['sourceGrid']} differs from source grid {args.source_grid}",
            flush=True,
        )

    bases = [int(piece) for piece in args.bases.replace(",", " ").split()]
    ratios = parse_floats(args.global_constraint_ratios)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    tail_budget = bridge_tail_budget(args.bridge_json)
    budget_fraction = mp.mpf(str(tail_budget["fraction"])) if tail_budget else None

    print("Full-theta source-inactive Schur-tail certificate", flush=True)
    print("  basis ratio inactive/top continuum/top active_trace active_ker full_ker", flush=True)
    rows = []
    for basis in bases:
        for ratio in ratios:
            row = compute_case(args, basis, ratio, e_derivs)
            rows.append(row)
            print(
                f"  {basis:5d} {ratio:5.3f} "
                f"{fmt(mp.mpf(row['finiteInactiveFracOfTop']), 8):>12} "
                f"{fmt(mp.mpf(row['continuumInactiveFracUpper']), 8):>13} "
                f"{fmt(mp.mpf(row['activeTraceMinNormed']), 8):>12} "
                f"{fmt(mp.mpf(row['activeTraceKernelSourceFrac']), 8):>10} "
                f"{fmt(mp.mpf(row['fullTraceKernelSourceFrac']), 8):>10}",
                flush=True,
            )

    worst_continuum = max(mp.mpf(row["continuumInactiveFracUpper"]) for row in rows)
    worst_active_kernel = max(mp.mpf(row["activeTraceKernelSourceFrac"]) for row in rows)
    worst_full_kernel = max(mp.mpf(row["fullTraceKernelSourceFrac"]) for row in rows)
    min_active_trace = min(mp.mpf(row["activeTraceMinNormed"]) for row in rows)
    budget_pass = bool(budget_fraction is not None and worst_continuum < budget_fraction)
    sampled_pass = bool(worst_active_kernel == 0 and min_active_trace > 0)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceGrid": args.source_grid,
        "fullNmax": args.full_nmax,
        "activeCount": args.active_count,
        "quadratureJson": args.quadrature_json,
        "continuumErrorBound": f(mp.mpf(args.continuum_error)),
        "bridgeJson": args.bridge_json,
        "finiteSchurTailBudgetDiagnostic": tail_budget,
        "worstContinuumInactiveFracUpper": f(worst_continuum),
        "worstActiveTraceKernelSourceFrac": f(worst_active_kernel),
        "worstFullTraceKernelSourceFrac": f(worst_full_kernel),
        "minActiveTraceEigenvalueNormed": f(min_active_trace),
        "sampledActiveKernelExclusionPasses": sampled_pass,
        "finiteBudgetDiagnosticPasses": budget_pass,
        "globalSchurTailTheoremClosed": False,
        "rows": rows,
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Finite full-Phi source-inactive Schur-tail certificate.  The "
            "continuum inactive bound uses lambda_3(S_{8,h}) plus the certified "
            "full-Phi source quadrature/tail error.  The budget comparison is "
            "a finite diagnostic; the continuum high-frequency Hardy/Schur "
            "tail theorem remains the analytic target."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  worst continuum inactive/top <= {fmt(worst_continuum, 12)}")
    if budget_fraction is not None:
        print(f"  finite budget diagnostic {fmt(budget_fraction, 12)} pass={budget_pass}")
    print(f"  sampled active kernel exclusion pass={sampled_pass}")
    print(f"  global Schur tail theorem closed=False")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
