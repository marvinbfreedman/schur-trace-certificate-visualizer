#!/usr/bin/env python3
r"""Interval/refinement enclosure for the endpoint Krawczyk coefficients.

The previous ``endpoint_coefficient_ball_enclosure.py`` could be run in a
``--rounding-only`` mode, which only propagated an assumed sample ball through
the finite Krawczyk algebra.  This script removes that shortcut for the finite
coefficient-input layer:

* rebuild the scaled companion matrices and active endpoint boundary rows at a
  ladder of stricter confluent settings;
* use the refinement differences, plus an explicit tail model, as an entrywise
  interval radius around the strictest computed center;
* compare those radii with the Krawczyk capacities from
  ``endpoint_riccati_krawczyk_collocation.json``;
* state the exact proof dependency cleanly.

This is a finite endpoint-flow coefficient certificate.  It does not pretend
that refinement alone proves the analytic quadrature tail for the confluent
integral.  When the geometric-tail check is closed, the remaining formal item
is exactly the analytic theorem that the chosen refinement ladder dominates the
true confluent quadrature/eigenrow error.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    companion_matrix,
    status,
)
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    active_endpoint_covectors,
    aligned_derivative_cache,
    nearest_derivs,
    node_key,
    scaling_diag,
)
from endpoint_coefficient_ball_enclosure import (  # noqa: E402
    compare_objects,
)
from endpoint_flow_chebyshev_center import cheb_lobatto_nodes  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_level(text: str) -> tuple[int, int]:
    if ":" in text:
        dps_text, order_text = text.split(":", 1)
        return int(dps_text), int(order_text)
    value = int(text)
    return value, value


def parse_levels(text: str) -> list[tuple[int, int]]:
    levels = [parse_level(piece) for piece in text.replace(",", " ").split()]
    if len(levels) < 2:
        raise SystemExit("--levels must contain at least two entries")
    return levels


def level_tag(idx: int, total: int, dps: int, matrix_order: int, has_base_cache: bool) -> str:
    if idx == 0 and has_base_cache:
        return "base"
    if idx == total - 1 and dps >= 90 and matrix_order >= 90:
        return "strict"
    return f"level_{dps}_{matrix_order}"


def args_variant(args, dps: int, matrix_order: int, cache_suffix: str):
    out = SimpleNamespace(**vars(args))
    out.dps = dps
    out.matrix_order = matrix_order
    if cache_suffix == "base" and args.base_cache_dir:
        out.cache_dir = args.base_cache_dir
    else:
        out.cache_dir = str(Path(args.cache_dir) / cache_suffix)
    return out


def max_abs_matrix(mats: list[mp.matrix]) -> mp.mpf:
    return max(
        abs(mat[i, j])
        for mat in mats
        for i in range(mat.rows)
        for j in range(mat.cols)
    )


def matrix_max_abs(mat: mp.matrix) -> mp.mpf:
    return max(abs(mat[i, j]) for i in range(mat.rows) for j in range(mat.cols))


def build_objects_fixed_active(args, dps: int, matrix_order: int, cache_suffix: str, reference=None):
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
    if reference is None:
        base, active_modes, active_idx = active_setup(local_args, e_s0)
        reference = {
            "polys": base["polys"],
            "activeModes": active_modes,
            "activeIndices": active_idx,
        }
    c_left = active_endpoint_covectors(
        reference["polys"],
        reference["activeModes"],
        left,
        e_left,
        args.jet_order,
    )
    c_right = active_endpoint_covectors(
        reference["polys"],
        reference["activeModes"],
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
        "activeIndices": reference["activeIndices"],
        "minGap": min_gap,
        "minCos": min_cos,
        "derivativeRows": derivative_rows,
        "fixedActiveReference": reference,
        "scaledCompanionMaxAbs": max(max_abs_matrix(left_companions), max_abs_matrix(right_companions)),
        "boundaryRowMaxAbs": max(matrix_max_abs(c_left), matrix_max_abs(c_right)),
    }


def as_float_dict(row: dict) -> dict:
    out = {}
    for key, value in row.items():
        out[key] = f(value) if isinstance(value, mp.mpf) else value
    return out


def ratio_or_none(numer: mp.mpf, denom: mp.mpf):
    if denom == 0:
        return None
    return numer / denom


def geometric_tail(last_diff: mp.mpf, previous_diff: mp.mpf, ratio_cap: mp.mpf):
    if previous_diff <= 0:
        return mp.inf, mp.inf, False
    observed = last_diff / previous_diff
    q = max(observed, ratio_cap)
    if q >= 1:
        return mp.inf, observed, False
    return last_diff * q / (1 - q), observed, True


def last_difference_radius(
    diffs: list[dict],
    key: str,
    safety: mp.mpf,
    tail_mode: str,
    tail_ratio_cap: mp.mpf,
) -> tuple[mp.mpf, mp.mpf, mp.mpf | None, bool]:
    last = mp.mpf(diffs[-1][key])
    if tail_mode == "none":
        return safety * last, mp.mpf("0"), None, False
    if len(diffs) < 2:
        return safety * last, mp.mpf("0"), None, False
    previous = mp.mpf(diffs[-2][key])
    tail, observed, ok = geometric_tail(last, previous, tail_ratio_cap)
    return safety * (last + tail), tail, observed, ok


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--krawczyk-json", default="endpoint_riccati_krawczyk_consequence_theorem.json")
    parser.add_argument("--levels", default="70:70 90:90")
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--state-scaling", choices=["raw", "taylor"], default="taylor")
    parser.add_argument("--safety-factor", default="16")
    parser.add_argument("--tail-mode", choices=["none", "geometric"], default="none")
    parser.add_argument("--tail-ratio-cap", default="1e-6")
    parser.add_argument("--pilot-nodes", type=int, default=0)
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
    parser.add_argument("--json-out", default="endpoint_coefficient_interval_enclosure.json")
    args = parser.parse_args()

    levels = parse_levels(args.levels)
    kraw = load_json(args.krawczyk_json)
    coeff_capacity = mp.mpf(str(kraw["coefficientRadiusCapacityScaledCompanion"]))
    boundary_capacity = mp.mpf(str(kraw["boundaryRadiusCapacityUnscaledRows"]))
    simultaneous_fraction = mp.mpf(str(kraw.get("simultaneousCapacityFraction", "0")))
    safety = mp.mpf(args.safety_factor)
    tail_ratio_cap = mp.mpf(args.tail_ratio_cap)

    print(
        "Endpoint coefficient interval/refinement enclosure "
        f"levels={levels} tail={args.tail_mode}",
        flush=True,
    )

    objects = []
    fixed_reference = None
    for idx, (dps, matrix_order) in enumerate(levels):
        tag = level_tag(idx, len(levels), dps, matrix_order, bool(args.base_cache_dir))
        level_objects = build_objects_fixed_active(
            args,
            dps,
            matrix_order,
            tag,
            fixed_reference,
        )
        if fixed_reference is None:
            fixed_reference = level_objects["fixedActiveReference"]
        objects.append(
            {
                "level": {"dps": dps, "matrixOrder": matrix_order, "tag": tag},
                "objects": level_objects,
            }
        )

    diffs = []
    for left, right in zip(objects[:-1], objects[1:]):
        diff = compare_objects(left["objects"], right["objects"])
        row = {
            "from": left["level"],
            "to": right["level"],
            **as_float_dict(diff),
        }
        diffs.append(row)

    coeff_radius, coeff_tail, coeff_observed_ratio, coeff_tail_ok = last_difference_radius(
        diffs,
        "scaledCompanionMaxDiff",
        safety,
        args.tail_mode,
        tail_ratio_cap,
    )
    boundary_radius, boundary_tail, boundary_observed_ratio, boundary_tail_ok = last_difference_radius(
        diffs,
        "boundaryRowMaxDiff",
        safety,
        args.tail_mode,
        tail_ratio_cap,
    )
    coeff_pass = coeff_radius < coeff_capacity
    boundary_pass = boundary_radius < boundary_capacity
    within_independent_caps = coeff_pass and boundary_pass
    simultaneous_pass = (
        simultaneous_fraction > 0
        and coeff_radius <= simultaneous_fraction * coeff_capacity
        and boundary_radius <= simultaneous_fraction * boundary_capacity
    )
    tail_model_closed = (
        args.tail_mode == "geometric"
        and coeff_tail_ok
        and boundary_tail_ok
        and coeff_observed_ratio is not None
        and boundary_observed_ratio is not None
    )
    finite_krawczyk_budget_closed = within_independent_caps and simultaneous_pass
    exact_finite_closed = finite_krawczyk_budget_closed and tail_model_closed

    strict_level = objects[-1]["level"]
    strict_objects = objects[-1]["objects"]
    data = {
        "theoremName": "endpoint coefficient interval enclosure",
        "sourceKrawczykJson": args.krawczyk_json,
        "order": args.order,
        "stateScaling": args.state_scaling,
        "constraintInterval": [f(mp.mpf(args.constraint_min)), f(mp.mpf(args.constraint_max))],
        "s0": f(mp.mpf(args.s0)),
        "levels": [item["level"] for item in objects],
        "strictestLevel": strict_level,
        "activeCoordinatePolicy": (
            "fixed from the first/base Krawczyk active source basis; strict "
            "boundary rows are projected into the same active coordinates "
            "before entrywise comparison"
        ),
        "safetyFactor": f(safety),
        "tailMode": args.tail_mode,
        "tailRatioCap": f(tail_ratio_cap),
        "activeIndices": strict_objects["activeIndices"],
        "strictTraceFieldMinGap": f(strict_objects["minGap"]),
        "strictTraceFieldMinConsecutiveCos": f(strict_objects["minCos"]),
        "refinementDiffs": diffs,
        "scaledCompanionIntervalRadius": f(coeff_radius),
        "scaledCompanionTailAllowance": f(coeff_tail),
        "scaledCompanionObservedTailRatio": f(coeff_observed_ratio) if coeff_observed_ratio is not None else None,
        "scaledCompanionCapacity": f(coeff_capacity),
        "scaledCompanionRadiusOverCapacity": f(coeff_radius / coeff_capacity),
        "boundaryRowIntervalRadius": f(boundary_radius),
        "boundaryRowTailAllowance": f(boundary_tail),
        "boundaryRowObservedTailRatio": f(boundary_observed_ratio) if boundary_observed_ratio is not None else None,
        "boundaryRowCapacity": f(boundary_capacity),
        "boundaryRowRadiusOverCapacity": f(boundary_radius / boundary_capacity),
        "simultaneousCapacityFraction": f(simultaneous_fraction),
        "simultaneousCoefficientFractionUsed": f(coeff_radius / coeff_capacity),
        "simultaneousBoundaryFractionUsed": f(boundary_radius / boundary_capacity),
        "scaledCompanionIntervalStatus": status(
            "scaled companion interval below Krawczyk capacity",
            bool(coeff_pass),
            (
                "Closed when the safety-inflated refinement interval radius "
                "for every scaled companion entry is smaller than the "
                "Krawczyk coefficient capacity."
            ),
        ),
        "boundaryRowIntervalStatus": status(
            "active boundary-row interval below Krawczyk capacity",
            bool(boundary_pass),
            (
                "Closed when the safety-inflated refinement interval radius "
                "for every active endpoint boundary-row entry is smaller than "
                "the Krawczyk boundary-row capacity."
            ),
        ),
        "simultaneousKrawczykInputStatus": status(
            "simultaneous Krawczyk input below half-chart budget",
            bool(simultaneous_pass),
            (
                "Closed when both coefficient and boundary intervals fit "
                "inside the simultaneous fraction reported by the finite "
                "Krawczyk certificate."
            ),
        ),
        "refinementTailModelStatus": status(
            "refinement tail model",
            bool(tail_model_closed),
            (
                "Closed only with tail-mode=geometric when the last two "
                "refinement differences support a convergent geometric tail. "
                "With tail-mode=none the Krawczyk budget comparison is still "
                "valid for the displayed finite refinement radius, but this "
                "status remains open as an analytic infinite-tail proof."
            ),
        ),
        "finiteKrawczykBudgetStatus": status(
            "finite Krawczyk budget comparison",
            bool(finite_krawczyk_budget_closed),
            (
                "Closed when the displayed finite/refinement radii for the "
                "fixed-active-coordinate coefficients fit inside both the "
                "independent and simultaneous Krawczyk budgets."
            ),
        ),
        "exactFiniteKrawczykCoefficientInputStatus": status(
            "exact finite Krawczyk coefficient input",
            bool(exact_finite_closed),
            (
                "Closed only after the finite Krawczyk budget comparison and "
                "the refinement tail model are both closed.  The remaining "
                "external proof obligation is to justify the refinement tail "
                "from analytic quadrature/eigenrow perturbation estimates."
            ),
        ),
        "fullRankShortcutStatus": status(
            "endpoint full-rank shortcut from coefficient enclosure",
            bool(exact_finite_closed and simultaneous_pass),
            (
                "If the coefficient intervals are exact, the Krawczyk theorem "
                "keeps the endpoint map in the persistent p_01 chart, so the "
                "active endpoint map has full row rank and the Green BVP has "
                "no active compatibility obstruction."
            ),
        ),
        "remainingFormalItem": (
            "Replace the refinement/tail model by a theorem bounding the "
            "confluent quadrature, r-tail, and eigenrow perturbation errors "
            "for every node used in the order-11 collocation system."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint coefficient interval enclosure")
    print(f"  scaled companion radius: {fmt(coeff_radius, 12)} / cap {fmt(coeff_capacity, 12)}")
    print(f"  boundary row radius: {fmt(boundary_radius, 12)} / cap {fmt(boundary_capacity, 12)}")
    print(f"  simultaneous fraction used: coeff={fmt(coeff_radius/coeff_capacity, 12)} boundary={fmt(boundary_radius/boundary_capacity, 12)}")
    print(f"  independent-cap pass: {within_independent_caps}")
    print(f"  simultaneous pass: {simultaneous_pass}")
    print(f"  exact finite input closed: {exact_finite_closed}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
