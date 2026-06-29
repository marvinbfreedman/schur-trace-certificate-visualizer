#!/usr/bin/env python3
r"""Magnus/Taylor-scaled center for the endpoint adjoint flow.

This is the next endpoint-flow center after Chebyshev collocation.  It solves

    Y'(s)=A(s)Y(s)

with a fourth-order two-point Gauss-Magnus step:

    Omega = h/2(A1+A2) + sqrt(3) h^2/12 [A2,A1].

Each step is evaluated in a local Taylor-scaled derivative basis before taking
``expm(Omega)`` and then converted back to raw derivative coordinates.  The
coefficient matrices A(s) are computed from exact confluent trace derivatives,
not finite differences of the moving eigenrow.

This is still a numerical center and refinement test, not an interval proof.
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
    endpoint_map,
    left_obstruction,
    mat_frobenius,
    matrix_rank_gram,
    min_norm_solve,
    row_norm,
    status,
)
from endpoint_flow_chebyshev_center import (  # noqa: E402
    compare_centers,
    exact_derivative_at,
    matrix_from_rows,
    source_endpoint_rows,
)
from global_trace_observability_gap import f, fmt  # noqa: E402


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def flip_derivs(e_derivs):
    return [[-value for value in row] for row in e_derivs]


def vec_dot(a, b) -> mp.mpf:
    return mp.fsum(a[i] * b[i] for i in range(len(a)))


def scaling_diag(radius: mp.mpf, dim: int) -> mp.matrix:
    scale = mp.matrix(dim, dim)
    r = max(abs(radius), mp.mpf("1e-30"))
    for j in range(dim):
        scale[j, j] = r**j / mp.factorial(j)
    return scale


def scaled_step_matrix(amat: mp.matrix, scale: mp.matrix) -> mp.matrix:
    return scale * amat * (scale ** -1)


def exact_e_derivs(args, s: mp.mpf, max_q: int):
    vals, e_derivs, lam_derivs = exact_derivative_at(args, s, max_q)
    gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
    return vals, e_derivs, lam_derivs, gap


def companion_at(args, s: mp.mpf) -> mp.matrix:
    _vals, e_derivs, _lam, _gap = exact_e_derivs(args, s, args.jet_order - 1)
    return companion_matrix(e_derivs, args.jet_order - 1)


def magnus_step(args, start: mp.mpf, end: mp.mpf) -> mp.matrix:
    dim = args.jet_order - 1
    h = end - start
    c = mp.sqrt(mp.mpf("3")) / 6
    s1 = start + (mp.mpf("0.5") - c) * h
    s2 = start + (mp.mpf("0.5") + c) * h
    scale = scaling_diag(abs(h), dim)
    scale_inv = scale ** -1
    a1 = scaled_step_matrix(companion_at(args, s1), scale)
    a2 = scaled_step_matrix(companion_at(args, s2), scale)
    comm = a2 * a1 - a1 * a2
    omega = (h / 2) * (a1 + a2) + (mp.sqrt(mp.mpf("3")) * h * h / 12) * comm
    return scale_inv * mp.expm(omega) * scale


def magnus_flow(args, start: mp.mpf, end: mp.mpf, steps: int) -> mp.matrix:
    dim = args.jet_order - 1
    flow = mp.eye(dim)
    h = (end - start) / steps
    current = start
    for _ in range(steps):
        nxt = current + h
        flow = magnus_step(args, current, nxt) * flow
        current = nxt
    return flow


def align_to(reference, e_derivs):
    if vec_dot(reference[0], e_derivs[0]) < 0:
        return flip_derivs(e_derivs)
    return e_derivs


def build_center(args, steps: int) -> dict:
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    full_q = max(args.max_trace_q, args.jet_order - 1)
    vals_s0, e_s0, _lam_s0, gap_s0 = exact_e_derivs(args, s0, full_q)
    _vals_left, e_left, _lam_left, gap_left = exact_e_derivs(args, left, args.jet_order - 1)
    _vals_right, e_right, _lam_right, gap_right = exact_e_derivs(args, right, args.jet_order - 1)
    e_left = align_to(e_s0, e_left)
    e_right = align_to(e_s0, e_right)

    flow_left = magnus_flow(args, s0, left, steps)
    flow_right = magnus_flow(args, s0, right, steps)

    base, active_modes, active_idx = active_setup(args, e_s0)
    endpoint_mat, _zero = endpoint_map(
        base["polys"],
        active_modes,
        e_left,
        e_right,
        flow_left,
        flow_right,
        left,
        right,
    )
    rank, eigs, rank_threshold = matrix_rank_gram(
        endpoint_mat * endpoint_mat.T,
        mp.mpf(args.trace_tol),
    )
    obstruction, obstruction_eigs = left_obstruction(endpoint_mat)
    obstruction_residual = endpoint_mat.T * obstruction
    obstruction_residual_norm = row_norm(
        [obstruction_residual[i] for i in range(obstruction_residual.rows)]
    )
    raw_rows = source_endpoint_rows(
        args,
        base,
        active_modes,
        e_left,
        e_right,
        flow_right,
        e_s0,
    )
    max_residual = mp.mpf("0")
    max_z_norm = mp.mpf("0")
    max_obstruction_pairing = mp.mpf("0")
    max_obstruction_relative_pairing = mp.mpf("0")
    rows = []
    for row in raw_rows:
        b_values = [mp.mpf(str(value)) for value in row["endpointVector"]]
        z, endpoint_residual, _gram = min_norm_solve(
            endpoint_mat,
            [-value for value in b_values],
        )
        residual_norm = row_norm([endpoint_residual[i] for i in range(endpoint_residual.rows)])
        z_norm = row_norm([z[i] for i in range(z.rows)])
        beta_before = row_norm(b_values)
        obstruction_pairing = mp.fsum(
            obstruction[i] * b_values[i] for i in range(len(b_values))
        )
        obstruction_relative = abs(obstruction_pairing) / max(mp.mpf("1"), beta_before)
        max_residual = max(max_residual, residual_norm)
        max_z_norm = max(max_z_norm, z_norm)
        max_obstruction_pairing = max(max_obstruction_pairing, abs(obstruction_pairing))
        max_obstruction_relative_pairing = max(
            max_obstruction_relative_pairing,
            obstruction_relative,
        )
        row.update(
            {
                "leftObstructionPairing": f(obstruction_pairing),
                "leftObstructionRelativePairing": f(obstruction_relative),
                "selectedZNorm": f(z_norm),
                "endpointResidualNorm": f(residual_norm),
                "endpointRelativeResidual": f(residual_norm / max(mp.mpf("1"), beta_before)),
                "selectedZ": [f(z[i]) for i in range(z.rows)],
            }
        )
        rows.append(row)

    sigmas = [mp.sqrt(max(mp.mpf("0"), value)) for value in eigs]
    return {
        "magnusStepsPerHalf": steps,
        "activeDim": len(active_idx),
        "activeIndices": active_idx,
        "traceFieldMinGap": f(min(gap_left, gap_s0, gap_right)),
        "endpointMapRank": rank,
        "endpointMapGramEigenvalues": [f(value) for value in eigs],
        "endpointMapRankThreshold": f(rank_threshold),
        "sigmaMin": f(min(sigmas)),
        "sigmaMax": f(max(sigmas)),
        "endpointMap": [
            [f(endpoint_mat[i, j]) for j in range(endpoint_mat.cols)]
            for i in range(endpoint_mat.rows)
        ],
        "leftObstructionVector": [f(obstruction[i]) for i in range(obstruction.rows)],
        "leftObstructionGramEigenvalues": [f(value) for value in obstruction_eigs],
        "leftObstructionResidualNorm": f(obstruction_residual_norm),
        "maxLeftObstructionPairing": f(max_obstruction_pairing),
        "maxLeftObstructionRelativePairing": f(max_obstruction_relative_pairing),
        "endpointMapFrobenius": f(mat_frobenius(endpoint_mat)),
        "flowLeftFrobenius": f(mat_frobenius(flow_left)),
        "flowRightFrobenius": f(mat_frobenius(flow_right)),
        "maxEndpointResidualNorm": f(max_residual),
        "maxSelectedZNorm": f(max_z_norm),
        "rows": rows,
        "statuses": {
            "activeEndpointMapFullRowRank": status(
                "active endpoint map full row rank",
                rank == len(active_idx),
                "Magnus endpoint center has full sampled active row rank.",
            ),
            "actualEndpointRhsInRange": status(
                "actual endpoint RHS in Range(M)",
                max_residual < mp.mpf("1e-20"),
                "For the actual source rows, endpoint vectors b(d) are in the computed endpoint-map range.",
            ),
            "activeEndpointKilled": status(
                "active endpoint concomitant killed",
                max_residual < mp.mpf("1e-20"),
                "Minimum-norm z solves Mz=-b(d) on the active source space.",
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", default="2 4 8")
    parser.add_argument("--rank-margin", default="2.77565955249e6")
    parser.add_argument("--cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
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
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
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
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="endpoint_flow_magnus_center.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    steps_list = parse_ints(args.steps)
    rank_margin = mp.mpf(args.rank_margin)
    centers = []
    for steps in steps_list:
        print(f"Building Magnus endpoint center steps={steps}", flush=True)
        center = build_center(args, steps)
        centers.append(center)
        print(
            f"  steps={steps} rank={center['endpointMapRank']} "
            f"||M||_F={fmt(mp.mpf(center['endpointMapFrobenius']), 12)} "
            f"sigma_min={fmt(mp.mpf(center['sigmaMin']), 12)} "
            f"endpoint_res={fmt(mp.mpf(center['maxEndpointResidualNorm']), 8)}",
            flush=True,
        )

    comparisons = []
    for left, right in zip(centers[:-1], centers[1:]):
        item = compare_centers(
            {
                "chebyshevOrder": left["magnusStepsPerHalf"],
                "endpointMap": left["endpointMap"],
            },
            {
                "chebyshevOrder": right["magnusStepsPerHalf"],
                "endpointMap": right["endpointMap"],
            },
            rank_margin,
        )
        item["leftSteps"] = item.pop("leftOrder")
        item["rightSteps"] = item.pop("rightOrder")
        comparisons.append(item)
    stable = bool(comparisons and all(row["passesRankMargin"] for row in comparisons))
    data = {
        "theoremName": "endpoint flow Magnus center",
        "interval": [f(mp.mpf(args.constraint_min)), f(mp.mpf(args.constraint_max))],
        "s0": f(mp.mpf(args.s0)),
        "steps": steps_list,
        "rankMargin": f(rank_margin),
        "centers": centers,
        "comparisons": comparisons,
        "magnusCenterStabilityStatus": status(
            "Magnus endpoint center stability",
            stable,
            (
                "Closed only if successive Magnus endpoint maps differ by "
                "less than the rank margin.  This is a numerical center test, "
                "not an interval proof."
            ),
        ),
        "nextTarget": (
            "If Magnus centers stabilize, wrap the stepper in interval/ball "
            "arithmetic.  If not, return to analytic obstruction annihilation."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint flow Magnus center")
    for row in comparisons:
        print(
            f"  {row['leftSteps']}->{row['rightSteps']} "
            f"diff={row['frobeniusDifference']:.12e} "
            f"diff/margin={row['differenceOverMargin']:.12e} "
            f"pass={row['passesRankMargin']}",
            flush=True,
        )
    print(f"  stable below rank margin: {stable}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
