#!/usr/bin/env python3
r"""Extract the sampled trace-to-source representation kernel.

The finite active-range identity

    E_active = C_sample R_global

is the Galerkin analogue of the desired continuum representation

    E_active f = int_I K(a)^* Lambda_a(f) da.

This script constructs ``C_sample`` on the representative source-active high
block and records each row of ``C_sample`` as a coefficient profile over the
sampled trace centers.  Bounded, stable profiles are the numerical footprint of
an L^1/L^2 trace-resolution kernel.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import constraints_for, local_case  # noqa: E402
from global_trace_active_range_inclusion import inverse_spd, op_norm  # noqa: E402
from global_trace_observability_gap import (  # noqa: E402
    child_args,
    f,
    fmt,
    frob_norm,
    parse_ints,
    source_nodes,
    stack_source_rows,
)
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402


def profile_stats(values, centers):
    vals = [mp.mpf(v) for v in values]
    if not vals:
        return {
            "l1Grid": 0.0,
            "l2Grid": 0.0,
            "maxAbs": 0.0,
            "totalVariation": 0.0,
            "signedSum": 0.0,
        }
    if len(vals) == 1:
        weights = [mp.mpf("1")]
    else:
        weights = []
        for i in range(len(vals)):
            if i == 0:
                weights.append(abs(centers[1] - centers[0]) / 2)
            elif i == len(vals) - 1:
                weights.append(abs(centers[-1] - centers[-2]) / 2)
            else:
                weights.append(abs(centers[i + 1] - centers[i - 1]) / 2)
    l1 = mp.fsum(weights[i] * abs(vals[i]) for i in range(len(vals)))
    l2 = mp.sqrt(mp.fsum(weights[i] * abs(vals[i]) ** 2 for i in range(len(vals))))
    variation = mp.fsum(abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1))
    signed_sum = mp.fsum(weights[i] * vals[i] for i in range(len(vals)))
    return {
        "l1Grid": f(l1),
        "l2Grid": f(l2),
        "maxAbs": f(max(abs(v) for v in vals)),
        "totalVariation": f(variation),
        "signedIntegralGrid": f(signed_sum),
    }


def active_indices(source_vals, source_top, active_tol):
    floor = active_tol * max(mp.mpf("1"), source_top)
    return [i for i, val in enumerate(source_vals) if val > floor], floor


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument(
        "--global-constraint-rule",
        choices=["fixed", "ratio-floor", "ratio-round", "offset"],
        default="ratio-floor",
    )
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-8")
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
    parser.add_argument("--json-out", default="trace_to_source_kernel_profile.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    active_tol = mp.mpf(args.active_tol)
    max_local_constraints = constraints_for(args, args.basis, args.global_constraint_ratio, local=True)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    base = local_case(args, args.basis, e_derivs)
    nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs)
    base["sourceOnScaled"] = source_rows * base["scaledModes"]

    global_constraints = constraints_for(args, args.basis, args.global_constraint_ratio, local=False)
    gargs = make_qargs(child_args(args, args.basis, global_constraints))
    centers, R_global = trace_matrix(base["polys"], gargs)

    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], active_tol)
    active_basis = columns(base["sourceVecs"], active_idx)
    source_active = base["sourceOnScaled"] * active_basis
    trace_active = R_global * base["scaledModes"] * active_basis

    gram = trace_active.T * trace_active
    gvals, gkeep, gplus, gtol = inverse_spd(gram, args.trace_tol)
    coeff = source_active * gplus * trace_active.T
    residual = source_active - coeff * trace_active
    source_norm = frob_norm(source_active)
    residual_norm = frob_norm(residual)

    row_profiles = []
    for source_i, node in enumerate(nodes):
        for component in range(2):
            row = 2 * source_i + component
            coeffs = [coeff[row, j] for j in range(coeff.cols)]
            row_profiles.append(
                {
                    "sourceNode": f(node),
                    "component": "boundary" if component == 0 else "adjointEval",
                    "coefficients": [f(x) for x in coeffs],
                    "stats": profile_stats(coeffs, centers),
                }
            )

    source_singular = []
    if active_idx:
        for local_col, idx in enumerate(active_idx):
            source_singular.append(
                {
                    "sourceEigenIndex": idx,
                    "sourceEigenvalue": f(base["sourceVals"][idx]),
                    "relativeToTop": f(base["sourceVals"][idx] / base["sourceTop"]),
                }
            )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "cutoff": args.cutoff,
        "localConstraints": base["localConstraints"],
        "globalConstraints": global_constraints,
        "constraintMin": f(mp.mpf(args.constraint_min)),
        "constraintMax": f(mp.mpf(args.constraint_max)),
        "traceCenters": [f(x) for x in centers],
        "sourceNodes": [f(x) for x in nodes],
        "activeTol": f(active_tol),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "sourceTop": f(base["sourceTop"]),
        "sourceSingularDirections": source_singular,
        "traceRankOnActive": len(gkeep),
        "traceGramEigenvalues": [f(x) for x in gvals],
        "traceGramTol": f(gtol),
        "coefficientFrobenius": f(frob_norm(coeff)),
        "coefficientOperator": f(op_norm(coeff)),
        "sourceActiveFrobenius": f(source_norm),
        "rangeResidualFrobenius": f(residual_norm),
        "rangeResidualRelative": f(residual_norm / source_norm if source_norm else mp.mpf("0")),
        "rowProfiles": row_profiles,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "interpretation": (
            "Sampled trace-to-source representation.  For every source row "
            "listed here and every vector in the source-active high block, "
            "E_row f equals the displayed weighted sum of sampled Lambda_a(f), "
            "up to the reported residual."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(
        f"Trace-to-source kernel profile basis={args.basis} "
        f"global_constraints={global_constraints} active_dim={len(active_idx)}",
        flush=True,
    )
    print(
        f"  trace_rank={len(gkeep)} rel_res={fmt(residual_norm / source_norm if source_norm else 0, 12)} "
        f"coeff_frob={fmt(frob_norm(coeff), 12)} coeff_op={fmt(op_norm(coeff), 12)}",
        flush=True,
    )
    print("  source component      max|K|        L1_grid       variation", flush=True)
    for row in row_profiles:
        stats = row["stats"]
        print(
            f"  {fmt(mp.mpf(row['sourceNode']), 4):>6} {row['component']:<11} "
            f"{fmt(mp.mpf(stats['maxAbs']), 8):>12} "
            f"{fmt(mp.mpf(stats['l1Grid']), 8):>12} "
            f"{fmt(mp.mpf(stats['totalVariation']), 8):>12}",
            flush=True,
        )
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
