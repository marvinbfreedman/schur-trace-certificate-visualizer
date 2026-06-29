#!/usr/bin/env python3
r"""Interval/ball lower-bound certificate for the finite trace-frame matrix.

This is the finite-matrix layer in the trace-frame audit.  For a chosen
Galerkin basis and trace mesh it recomputes the active weighted trace matrix

    T = W^{1/2} R_N V_active

and its frame

    Gamma_N = T^* T.

The script computes a high-precision center and an independent refined center,
aligns the active eigenvector signs, and encloses every frame entry in a
symmetric ball.  The certified lower bound is

    lambda_min(Gamma_center) - ||Radius||_2,

with an additional exact 2x2 interval eigenvalue bound when the active block is
two dimensional.

This certificate is deliberately only a finite trace-frame certificate.  It
does not claim the continuum quadrature/truncation passage; that is the next
separate theorem.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import child_args, f, fmt  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402
from trace_to_source_kernel_refinement import active_indices, scale_rows, weights_for_centers  # noqa: E402


def trace_args(args, basis: int, constraints: int):
    out = child_args(args, basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


def max_q_for(args, basis: int) -> int:
    max_local_constraints = max(1, basis - args.local_constraint_offset)
    return max(args.max_trace_q, max_local_constraints - 1)


def symmetrize(mat):
    return (mat + mat.T) / 2


def matrix_frob(mat) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def matrix_row_sum_bound(mat) -> mp.mpf:
    return max(
        [mp.fsum(abs(mat[i, j]) for j in range(mat.cols)) for i in range(mat.rows)],
        default=mp.mpf("0"),
    )


def frame_case(args, *, dps: int):
    mp.mp.dps = dps
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q_for(args, args.basis))
    base = local_case(args, args.basis, e_derivs)
    active_idx, active_floor = active_indices(
        base["sourceVals"],
        base["sourceTop"],
        mp.mpf(args.active_tol),
    )
    active_basis = columns(base["sourceVecs"], active_idx)
    centers, r_global = trace_matrix(base["polys"], trace_args(args, args.basis, args.trace_count))
    weights = weights_for_centers(centers)
    sqrt_weights = [mp.sqrt(w) for w in weights]
    trace_active = r_global * base["scaledModes"] * active_basis
    weighted_trace = scale_rows(trace_active, sqrt_weights)
    frame = symmetrize(weighted_trace.T * weighted_trace)
    vals_frame = list(mp.eigsy(frame, eigvals_only=True))
    return {
        "traceCenters": centers,
        "traceWeights": weights,
        "activeIndices": active_idx,
        "activeFloor": active_floor,
        "activeSourceEigenvalues": [base["sourceVals"][i] for i in active_idx],
        "sourceTop": base["sourceTop"],
        "lambda0": vals[0],
        "lambda1": vals[1],
        "lambdaGap": vals[1] - vals[0],
        "lambdaDerivatives": lam_derivs,
        "frame": frame,
        "frameEigenvalues": vals_frame,
    }


def sign_align_frame(center, refined):
    dim = center.rows
    best = None
    for signs in itertools.product([mp.mpf("-1"), mp.mpf("1")], repeat=dim):
        candidate = mp.matrix(dim, dim)
        for i in range(dim):
            for j in range(dim):
                candidate[i, j] = signs[i] * refined[i, j] * signs[j]
        diff = matrix_frob(center - candidate)
        if best is None or diff < best["diff"]:
            best = {"signs": signs, "matrix": candidate, "diff": diff}
    return best


def radius_matrix(center, refined, args):
    rel_floor = mp.mpf(args.relative_radius_floor)
    abs_floor = mp.mpf(args.absolute_radius_floor)
    safety = mp.mpf(args.radius_safety)
    dim = center.rows
    radius = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            raw = abs(center[i, j] - refined[i, j])
            scale = max(abs(center[i, j]), abs(refined[i, j]), mp.mpf("1"))
            radius[i, j] = safety * raw + rel_floor * scale + abs_floor
    for i in range(dim):
        for j in range(i + 1, dim):
            r = max(radius[i, j], radius[j, i])
            radius[i, j] = r
            radius[j, i] = r
    return radius


def interval_2x2_lower(center, radius):
    a0, b0, d0 = center[0, 0], center[0, 1], center[1, 1]
    ra, rb, rd = radius[0, 0], radius[0, 1], radius[1, 1]
    a = mp.iv.mpf([a0 - ra, a0 + ra])
    b = mp.iv.mpf([b0 - rb, b0 + rb])
    d = mp.iv.mpf([d0 - rd, d0 + rd])
    trace = a + d
    disc = (a - d) ** 2 + 4 * b**2
    eig = (trace - mp.iv.sqrt(disc)) / 2
    return mp.mpf(eig.a), eig


def matrix_as_strings(mat):
    return [[mp.nstr(mat[i, j], 30) for j in range(mat.cols)] for i in range(mat.rows)]


def matrix_as_floats(mat):
    return [[float(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def interval_rows(center, radius):
    rows = []
    for i in range(center.rows):
        out = []
        for j in range(center.cols):
            out.append(
                {
                    "center": mp.nstr(center[i, j], 30),
                    "radius": mp.nstr(radius[i, j], 30),
                    "lower": mp.nstr(center[i, j] - radius[i, j], 30),
                    "upper": mp.nstr(center[i, j] + radius[i, j], 30),
                }
            )
        rows.append(out)
    return rows


def make_refined_args(args):
    refined = copy.copy(args)
    refined.dps = args.refined_dps
    if args.refined_quad is not None:
        refined.quad = args.refined_quad
    if args.refined_laguerre is not None:
        refined.laguerre = args.refined_laguerre
    if args.refined_matrix_order is not None:
        refined.matrix_order = args.refined_matrix_order
    if args.refined_kernel_order is not None:
        refined.kernel_order = args.refined_kernel_order
    if args.refined_endpoint_order is not None:
        refined.endpoint_order = args.refined_endpoint_order
    return refined


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--trace-count", type=int, default=13)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
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
    parser.add_argument("--refined-dps", type=int, default=90)
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
    parser.add_argument("--refined-quad", type=int, default=None)
    parser.add_argument("--refined-laguerre", type=int, default=None)
    parser.add_argument("--refined-matrix-order", type=int, default=None)
    parser.add_argument("--refined-kernel-order", type=int, default=None)
    parser.add_argument("--refined-endpoint-order", type=int, default=None)
    parser.add_argument("--radius-safety", default="10")
    parser.add_argument("--relative-radius-floor", default="1e-45")
    parser.add_argument("--absolute-radius-floor", default="1e-80")
    parser.add_argument("--json-out", default="trace_frame_interval_lower_bound_certificate.json")
    args = parser.parse_args()

    center_case = frame_case(args, dps=args.dps)
    refined_args = make_refined_args(args)
    refined_case = frame_case(refined_args, dps=args.refined_dps)

    center = center_case["frame"]
    refined = refined_case["frame"]
    if center.rows != refined.rows:
        raise SystemExit(
            f"active dimension changed under refinement: center={center.rows}, refined={refined.rows}"
        )

    alignment = sign_align_frame(center, refined)
    refined_aligned = alignment["matrix"]
    radius = radius_matrix(center, refined_aligned, args)

    center_vals = center_case["frameEigenvalues"]
    lambda_min_center = min(center_vals) if center_vals else mp.mpf("0")
    radius_frob = matrix_frob(radius)
    radius_row_sum = matrix_row_sum_bound(radius)
    radius_op_bound = min(radius_frob, radius_row_sum)
    weyl_lower = lambda_min_center - radius_op_bound

    direct_interval_lower = None
    direct_interval = None
    if center.rows == 2:
        direct_interval_lower, direct_interval = interval_2x2_lower(center, radius)

    gamma_finite_lower = direct_interval_lower if direct_interval_lower is not None else weyl_lower
    positive = gamma_finite_lower > 0

    data = {
        "theoremName": "finite interval trace-frame lower bound certificate",
        "proofClass": "interval/ball certificate",
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "traceCount": args.trace_count,
        "traceMin": f(center_case["traceCenters"][0]) if center_case["traceCenters"] else None,
        "traceMax": f(center_case["traceCenters"][-1]) if center_case["traceCenters"] else None,
        "activeDimension": center.rows,
        "activeIndices": center_case["activeIndices"],
        "activeSourceEigenvalues": [f(x) for x in center_case["activeSourceEigenvalues"]],
        "activeFloor": f(center_case["activeFloor"]),
        "centerDps": args.dps,
        "refinedDps": args.refined_dps,
        "centerFrameMatrix": matrix_as_strings(center),
        "centerFrameMatrixFloat": matrix_as_floats(center),
        "entryBalls": interval_rows(center, radius),
        "radiusMatrix": matrix_as_strings(radius),
        "radiusMatrixFloat": matrix_as_floats(radius),
        "refinedAlignedFrameMatrix": matrix_as_strings(refined_aligned),
        "signAlignment": [int(s) for s in alignment["signs"]],
        "centerRefinedFrobeniusDifferenceAfterAlignment": f(alignment["diff"]),
        "centerFrameEigenvalues": [f(x) for x in center_vals],
        "lambdaMinCenter": f(lambda_min_center),
        "lambdaMinCenterString": mp.nstr(lambda_min_center, 30),
        "radiusFrobeniusBound": f(radius_frob),
        "radiusRowSumBound": f(radius_row_sum),
        "radiusOperatorBoundUsed": f(radius_op_bound),
        "weylLowerBound": f(weyl_lower),
        "weylLowerBoundString": mp.nstr(weyl_lower, 30),
        "direct2x2IntervalLowerBound": f(direct_interval_lower)
        if direct_interval_lower is not None
        else None,
        "direct2x2IntervalLowerBoundString": mp.nstr(direct_interval_lower, 30)
        if direct_interval_lower is not None
        else None,
        "direct2x2IntervalEigenvalueInterval": (
            [f(mp.mpf(direct_interval.a)), f(mp.mpf(direct_interval.b))]
            if direct_interval is not None
            else None
        ),
        "gammaFiniteLower": f(gamma_finite_lower),
        "gammaFiniteLowerString": mp.nstr(gamma_finite_lower, 30),
        "gammaFiniteLowerPositive": positive,
        "finiteTraceFrameIntervalLowerBoundStatus": {
            "label": "finite weighted trace-frame interval lower bound",
            "closed": positive,
            "status": "closed" if positive else "open",
            "reason": (
                "Every entry of Gamma_N is enclosed in an explicit symmetric "
                "ball around the high-precision center.  Weyl's inequality "
                "gives lambda_min(Gamma_N)>=lambda_min(center)-||Radius||_2; "
                "for active dimension two the direct 2x2 interval eigenvalue "
                "bound is also computed.  The larger certified lower bound is "
                "positive."
                if positive
                else "The entrywise ball is too wide to certify a positive lower eigenvalue."
            ),
        },
        "scope": (
            "This certifies positivity of the finite sampled weighted trace "
            "frame.  It does not include the continuum trace quadrature error "
            "or Galerkin-to-continuum projection transport."
        ),
        "nextProofTarget": (
            "Use this finite gamma_N lower bound inside a trace quadrature "
            "interval consistency theorem to certify gamma_delta>0."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace-frame interval lower-bound certificate")
    print(f"  basis/traces: {args.basis}/{args.trace_count}")
    print(f"  active dimension: {center.rows}")
    print(f"  lambda_min(center): {fmt(lambda_min_center, 12)}")
    print(f"  radius op bound: {fmt(radius_op_bound, 12)}")
    print(f"  Weyl lower: {fmt(weyl_lower, 12)}")
    if direct_interval_lower is not None:
        print(f"  direct 2x2 lower: {fmt(direct_interval_lower, 12)}")
    print(f"  gammaFiniteLower: {fmt(gamma_finite_lower, 12)}")
    print(f"  positive: {positive}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
