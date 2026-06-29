#!/usr/bin/env python3
r"""Coefficient/boundary-row ball enclosure for the Krawczyk endpoint proof.

The Krawczyk reduction reports two concrete radii:

* an entrywise radius for the scaled companion matrices used in the row-flow
  collocation system;
* an entrywise radius for the active endpoint boundary rows.

This script recomputes those finite objects with a stricter confluent
integration setting and turns the refinement difference into an explicit ball
radius.  It is a certified finite-enclosure step conditional on the usual
validated-quadrature interpretation of the refinement ball; it does not yet
derive a symbolic quadrature tail bound for the confluent integral.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import active_setup, companion_matrix, status  # noqa: E402
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    active_endpoint_covectors,
    aligned_derivative_cache,
    nearest_derivs,
    node_key,
    scaling_diag,
)
from endpoint_flow_chebyshev_center import cheb_lobatto_nodes  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def max_abs_diff_matrix(a: mp.matrix, b: mp.matrix) -> mp.mpf:
    return max(abs(a[i, j] - b[i, j]) for i in range(a.rows) for j in range(a.cols))


def matrix_max_abs(a: mp.matrix) -> mp.mpf:
    return max(abs(a[i, j]) for i in range(a.rows) for j in range(a.cols))


def args_variant(args, dps: int, matrix_order: int, cache_suffix: str):
    out = SimpleNamespace(**vars(args))
    out.dps = dps
    out.matrix_order = matrix_order
    if cache_suffix == "base" and args.base_cache_dir:
        out.cache_dir = args.base_cache_dir
    else:
        out.cache_dir = str(Path(args.cache_dir) / cache_suffix)
    return out


def build_objects(args, dps: int, matrix_order: int, cache_suffix: str):
    local_args = args_variant(args, dps, matrix_order, cache_suffix)
    mp.mp.dps = dps
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    dim = args.jet_order - 1
    max_q = max(args.max_trace_q, dim)
    left_nodes = cheb_lobatto_nodes(args.order, left, s0)
    right_nodes = cheb_lobatto_nodes(args.order, s0, right)
    if args.pilot_nodes:
        keep = max(2, args.pilot_nodes)
        left_nodes = left_nodes[: keep - 1] + [left_nodes[-1]]
        right_nodes = [right_nodes[0]] + right_nodes[-(keep - 1) :]
    all_nodes = left_nodes + right_nodes
    print(
        f"  {cache_suffix}: dps={dps} matrix_order={matrix_order} "
        f"nodes={len(sorted({node_key(node): node for node in all_nodes}))}",
        flush=True,
    )
    cache, derivative_rows, min_gap, min_cos = aligned_derivative_cache(
        local_args,
        all_nodes,
        max_q,
    )
    scale_left = scaling_diag(left_nodes, args.state_scaling, dim)
    scale_left_inv = scale_left ** -1
    scale_right = scaling_diag(right_nodes, args.state_scaling, dim)
    scale_right_inv = scale_right ** -1
    left_companions = [
        scale_left * companion_matrix(cache[node_key(node)], dim) * scale_left_inv
        for node in left_nodes
    ]
    right_companions = [
        scale_right * companion_matrix(cache[node_key(node)], dim) * scale_right_inv
        for node in right_nodes
    ]

    e_s0 = nearest_derivs(cache, s0)
    e_left = nearest_derivs(cache, left)
    e_right = nearest_derivs(cache, right)
    base, active_modes, active_idx = active_setup(local_args, e_s0)
    c_left = active_endpoint_covectors(
        base["polys"],
        active_modes,
        left,
        e_left,
        args.jet_order,
    )
    c_right = active_endpoint_covectors(
        base["polys"],
        active_modes,
        right,
        e_right,
        args.jet_order,
    )
    return {
        "leftNodes": left_nodes,
        "rightNodes": right_nodes,
        "leftCompanions": left_companions,
        "rightCompanions": right_companions,
        "leftBoundaryRows": c_left,
        "rightBoundaryRows": c_right,
        "activeIndices": active_idx,
        "minGap": min_gap,
        "minCos": min_cos,
        "derivativeRows": derivative_rows,
    }


def compare_objects(base, strict):
    left_comp_diff = max(
        max_abs_diff_matrix(a, b)
        for a, b in zip(base["leftCompanions"], strict["leftCompanions"])
    )
    right_comp_diff = max(
        max_abs_diff_matrix(a, b)
        for a, b in zip(base["rightCompanions"], strict["rightCompanions"])
    )
    coeff_diff = max(left_comp_diff, right_comp_diff)
    left_boundary_diff = max_abs_diff_matrix(
        base["leftBoundaryRows"],
        strict["leftBoundaryRows"],
    )
    right_boundary_diff = max_abs_diff_matrix(
        base["rightBoundaryRows"],
        strict["rightBoundaryRows"],
    )
    boundary_diff = max(left_boundary_diff, right_boundary_diff)
    return {
        "leftCompanionMaxDiff": left_comp_diff,
        "rightCompanionMaxDiff": right_comp_diff,
        "scaledCompanionMaxDiff": coeff_diff,
        "leftBoundaryRowMaxDiff": left_boundary_diff,
        "rightBoundaryRowMaxDiff": right_boundary_diff,
        "boundaryRowMaxDiff": boundary_diff,
        "scaledCompanionMaxAbs": max(
            max(matrix_max_abs(mat) for mat in base["leftCompanions"]),
            max(matrix_max_abs(mat) for mat in base["rightCompanions"]),
        ),
        "boundaryRowMaxAbs": max(
            matrix_max_abs(base["leftBoundaryRows"]),
            matrix_max_abs(base["rightBoundaryRows"]),
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--krawczyk-json", default="endpoint_riccati_krawczyk_collocation.json")
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--base-dps", type=int, default=70)
    parser.add_argument("--strict-dps", type=int, default=90)
    parser.add_argument("--base-matrix-order", type=int, default=70)
    parser.add_argument("--strict-matrix-order", type=int, default=90)
    parser.add_argument("--safety-factor", default="16")
    parser.add_argument("--pilot-nodes", type=int, default=0)
    parser.add_argument("--rounding-only", action="store_true")
    parser.add_argument("--sample-rel-radius", default="1e-40")
    parser.add_argument("--sample-abs-radius", default="1e-45")
    parser.add_argument("--cache-dir", default=".endpoint_coeff_ball_cache")
    parser.add_argument("--base-cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=11)
    parser.add_argument("--active-tol", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
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
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="endpoint_coefficient_ball_enclosure.json")
    args = parser.parse_args()

    kraw = load_json(args.krawczyk_json)
    coeff_capacity = mp.mpf(str(kraw["coefficientRadiusCapacityScaledCompanion"]))
    boundary_capacity = mp.mpf(str(kraw["boundaryRadiusCapacityUnscaledRows"]))
    safety = mp.mpf(args.safety_factor)

    print(
        "Building endpoint coefficient/boundary ball enclosure "
        f"base=({args.base_dps},{args.base_matrix_order}) "
        f"strict=({args.strict_dps},{args.strict_matrix_order})",
        flush=True,
    )
    base = build_objects(args, args.base_dps, args.base_matrix_order, "base")
    if args.rounding_only:
        strict = None
        coeff_max_abs = max(
            max(matrix_max_abs(mat) for mat in base["leftCompanions"]),
            max(matrix_max_abs(mat) for mat in base["rightCompanions"]),
        )
        boundary_max_abs = max(
            matrix_max_abs(base["leftBoundaryRows"]),
            matrix_max_abs(base["rightBoundaryRows"]),
        )
        sample_rel = mp.mpf(args.sample_rel_radius)
        sample_abs = mp.mpf(args.sample_abs_radius)
        diff = {
            "leftCompanionMaxDiff": mp.mpf("0"),
            "rightCompanionMaxDiff": mp.mpf("0"),
            "scaledCompanionMaxDiff": max(sample_abs, sample_rel * coeff_max_abs),
            "leftBoundaryRowMaxDiff": mp.mpf("0"),
            "rightBoundaryRowMaxDiff": mp.mpf("0"),
            "boundaryRowMaxDiff": max(sample_abs, sample_rel * boundary_max_abs),
            "scaledCompanionMaxAbs": coeff_max_abs,
            "boundaryRowMaxAbs": boundary_max_abs,
        }
    else:
        strict = build_objects(args, args.strict_dps, args.strict_matrix_order, "strict")
        diff = compare_objects(base, strict)
    coeff_radius = safety * diff["scaledCompanionMaxDiff"]
    boundary_radius = safety * diff["boundaryRowMaxDiff"]
    coeff_pass = coeff_radius < coeff_capacity
    boundary_pass = boundary_radius < boundary_capacity
    data = {
        "theoremName": "endpoint coefficient and boundary-row ball enclosure",
        "sourceKrawczykJson": args.krawczyk_json,
        "order": args.order,
        "stateScaling": args.state_scaling,
        "baseDps": args.base_dps,
        "strictDps": args.strict_dps,
        "baseMatrixOrder": args.base_matrix_order,
        "strictMatrixOrder": args.strict_matrix_order,
        "matrixRmax": f(mp.mpf(args.matrix_rmax)),
        "safetyFactor": f(safety),
        "roundingOnly": bool(args.rounding_only),
        "sampleRelativeRadius": f(mp.mpf(args.sample_rel_radius)),
        "sampleAbsoluteRadius": f(mp.mpf(args.sample_abs_radius)),
        "activeIndices": base["activeIndices"],
        "baseTraceFieldMinGap": f(base["minGap"]),
        "strictTraceFieldMinGap": f(strict["minGap"]) if strict else None,
        "baseTraceFieldMinConsecutiveCos": f(base["minCos"]),
        "strictTraceFieldMinConsecutiveCos": f(strict["minCos"]) if strict else None,
        "scaledCompanionMaxDiff": f(diff["scaledCompanionMaxDiff"]),
        "leftCompanionMaxDiff": f(diff["leftCompanionMaxDiff"]),
        "rightCompanionMaxDiff": f(diff["rightCompanionMaxDiff"]),
        "boundaryRowMaxDiff": f(diff["boundaryRowMaxDiff"]),
        "leftBoundaryRowMaxDiff": f(diff["leftBoundaryRowMaxDiff"]),
        "rightBoundaryRowMaxDiff": f(diff["rightBoundaryRowMaxDiff"]),
        "scaledCompanionBallRadius": f(coeff_radius),
        "boundaryRowBallRadius": f(boundary_radius),
        "scaledCompanionCapacity": f(coeff_capacity),
        "boundaryRowCapacity": f(boundary_capacity),
        "scaledCompanionRadiusOverCapacity": f(coeff_radius / coeff_capacity),
        "boundaryRowRadiusOverCapacity": f(boundary_radius / boundary_capacity),
        "scaledCompanionMaxAbs": f(diff["scaledCompanionMaxAbs"]),
        "boundaryRowMaxAbs": f(diff["boundaryRowMaxAbs"]),
        "scaledCompanionEnclosureStatus": status(
            "scaled companion matrix ball below Krawczyk capacity",
            bool(coeff_pass),
            (
                "Closed when the safety-inflated refinement ball for every "
                "scaled companion entry is below the Krawczyk coefficient "
                "radius capacity."
            ),
        ),
        "boundaryRowEnclosureStatus": status(
            "active endpoint boundary-row ball below Krawczyk capacity",
            bool(boundary_pass),
            (
                "Closed when the safety-inflated refinement ball for every "
                "active endpoint boundary row entry is below the Krawczyk "
                "boundary radius capacity."
            ),
        ),
        "validatedKrawczykCoefficientInputStatus": status(
            "validated Krawczyk coefficient input",
            bool(coeff_pass and boundary_pass),
            (
                "This closes the finite coefficient-input layer conditional on "
                "the stated sample/refinement ball enclosure model.  A fully "
                "formal proof would replace this model with interval quadrature "
                "tail bounds for the confluent integral."
            ),
        ),
        "conclusion": (
            "The exact objects needed by the Krawczyk proof are now compared "
            "against stricter confluent computations.  Passing radii mean the "
            "scaled companion matrices and active endpoint rows are inside the "
            "Krawczyk budgets."
        ),
        "nextTarget": (
            "If this passes, state the endpoint full-rank theorem conditional "
            "on the validated confluent quadrature enclosure; if it fails, "
            "increase dps/matrix order or derive analytic quadrature tails."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint coefficient and boundary-row ball enclosure")
    print(f"  scaled companion diff: {fmt(diff['scaledCompanionMaxDiff'], 12)}")
    print(f"  scaled companion radius: {fmt(coeff_radius, 12)} / cap {fmt(coeff_capacity, 12)}")
    print(f"  boundary row diff: {fmt(diff['boundaryRowMaxDiff'], 12)}")
    print(f"  boundary row radius: {fmt(boundary_radius, 12)} / cap {fmt(boundary_capacity, 12)}")
    print(f"  companion pass: {coeff_pass}")
    print(f"  boundary pass: {boundary_pass}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
