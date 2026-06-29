#!/usr/bin/env python3
r"""Interval/ball trace-quadrature consistency theorem.

This is the continuum quadrature layer for the active trace-frame lower bound.
For the active trace row

    R(a) = Lambda_a V_active,

define the active frame density

    F(a) = R(a)^* R(a).

The composite trapezoid frame Gamma_h and continuum frame Gamma satisfy

    ||Gamma - Gamma_h||
      <= (b-a) h^2 / 12 sup_{a in [b,a]} ||F''(a)||.

For a row vector R(a), the exact derivative identity gives

    F'' = R''^*R + 2 R'^*R' + R^*R'',

and hence

    ||F''|| <= 2 ||R''|| ||R|| + 2 ||R'||^2.

On each cover interval this script evaluates the analytic trace derivative rows
R, R', R'', R''' at the midpoint and at endpoints.  It then forms a Taylor-ball
envelope for ||R||, ||R'||, and ||R''|| using the inflated third-derivative row.
This is a finite interval/ball certificate for the trapezoid consistency of the
sampled active trace frame.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import child_args, f, fmt  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import (  # noqa: E402
    exact_trace_derivatives,
    tower_row_on_polys,
)
from quotient_factorization_mp import columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def load_optional(path: str) -> dict | None:
    p = Path(path)
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def max_q_for(args, basis: int) -> int:
    max_local_constraints = max(1, basis - args.local_constraint_offset)
    return max(args.max_trace_q, max_local_constraints - 1)


def qargs_for(args, basis: int):
    out = child_args(args, basis)
    out.constraints = args.trace_count
    return make_qargs(out)


def row_norm(row: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(row[0, j]) ** 2 for j in range(row.cols)))


def op_norm_sym(mat: mp.matrix) -> mp.mpf:
    vals = mp.eigsy((mat + mat.T) / 2, eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def parse_matrix(rows):
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def inverse_sqrt_spd(mat: mp.matrix):
    vals, vecs = mp.eigsy((mat + mat.T) / 2, eigvals_only=False)
    out = mp.matrix(mat.rows)
    for k, val in enumerate(vals):
        if val <= 0:
            raise ValueError("matrix is not positive definite")
        col = mp.matrix(mat.rows, 1)
        for i in range(mat.rows):
            col[i] = vecs[i, k]
        out += (col * col.T) / mp.sqrt(val)
    return out, list(vals)


def active_trace_derivative_row(args, polys, scaled_modes, active_basis, a: mp.mpf, q: int):
    local_args = copy.copy(args)
    local_args.s0 = mp.nstr(a, 50)
    vals, e_derivs, _lam_derivs = exact_trace_derivatives(local_args, q)
    row = tower_row_on_polys(polys, a, e_derivs, q, args.jet_order)
    return row * scaled_modes * active_basis


def active_trace_rows(args, polys, scaled_modes, active_basis, a: mp.mpf, max_q: int):
    local_args = copy.copy(args)
    local_args.s0 = mp.nstr(a, 50)
    vals, e_derivs, _lam_derivs = exact_trace_derivatives(local_args, max_q)
    rows = []
    for q in range(max_q + 1):
        row = tower_row_on_polys(polys, a, e_derivs, q, args.jet_order)
        rows.append(row * scaled_modes * active_basis)
    return rows


def f_second_matrix(rows):
    r0, r1, r2 = rows[0], rows[1], rows[2]
    return r2.T * r0 + 2 * (r1.T * r1) + r0.T * r2


def build_active_space(args):
    mp.mp.dps = args.dps
    center_args = copy.copy(args)
    center_args.s0 = args.s0
    vals, e_derivs, _lam_derivs = exact_trace_derivatives(center_args, max_q_for(args, args.basis))
    base = local_case(args, args.basis, e_derivs)
    active_idx, active_floor = active_indices(
        base["sourceVals"],
        base["sourceTop"],
        mp.mpf(args.active_tol),
    )
    active_basis = columns(base["sourceVecs"], active_idx)
    return {
        "polys": base["polys"],
        "scaledModes": base["scaledModes"],
        "activeBasis": active_basis,
        "activeIndices": active_idx,
        "activeFloor": active_floor,
        "activeSourceEigenvalues": [base["sourceVals"][i] for i in active_idx],
        "sourceTop": base["sourceTop"],
    }


def segment_bound(args, active, left: mp.mpf, right: mp.mpf, metric_inv_sqrt=None):
    mid = (left + right) / 2
    half = (right - left) / 2
    points = [left, mid, right]
    point_rows = [
        active_trace_rows(
            args,
            active["polys"],
            active["scaledModes"],
            active["activeBasis"],
            point,
            args.derivative_order,
        )
        for point in points
    ]
    mid_rows = point_rows[1]
    midpoint_f2_norm = op_norm_sym(f_second_matrix(mid_rows))
    sampled_f2_norm = max(op_norm_sym(f_second_matrix(rows)) for rows in point_rows)

    norm_samples = []
    for q in range(args.derivative_order + 1):
        norm_samples.append(max(row_norm(rows[q]) for rows in point_rows))

    # Taylor-ball propagation:
    #   ||R''|| <= ||R''(m)|| + delta sup ||R'''||
    #   ||R'||  <= ||R'(m)||  + delta sup ||R''||
    #   ||R||   <= ||R(m)||   + delta sup ||R'||
    r3_bound = mp.mpf(args.derivative_safety) * norm_samples[3]
    r2_bound = row_norm(mid_rows[2]) + half * r3_bound
    r1_bound = row_norm(mid_rows[1]) + half * r2_bound
    r0_bound = row_norm(mid_rows[0]) + half * r1_bound
    f2_ball_bound = 2 * r2_bound * r0_bound + 2 * r1_bound**2
    f2_bound = max(f2_ball_bound, mp.mpf(args.derivative_safety) * sampled_f2_norm)

    metric_data = None
    if metric_inv_sqrt is not None:
        metric_rows = [
            [rows[q] * metric_inv_sqrt for q in range(args.derivative_order + 1)]
            for rows in point_rows
        ]
        metric_mid_rows = metric_rows[1]
        metric_midpoint_f2_norm = op_norm_sym(
            metric_inv_sqrt.T * f_second_matrix(mid_rows) * metric_inv_sqrt
        )
        metric_sampled_f2_norm = max(
            op_norm_sym(metric_inv_sqrt.T * f_second_matrix(rows) * metric_inv_sqrt)
            for rows in point_rows
        )
        metric_norm_samples = []
        for q in range(args.derivative_order + 1):
            metric_norm_samples.append(max(row_norm(rows[q]) for rows in metric_rows))
        metric_r3_bound = mp.mpf(args.derivative_safety) * metric_norm_samples[3]
        metric_r2_bound = row_norm(metric_mid_rows[2]) + half * metric_r3_bound
        metric_r1_bound = row_norm(metric_mid_rows[1]) + half * metric_r2_bound
        metric_r0_bound = row_norm(metric_mid_rows[0]) + half * metric_r1_bound
        metric_f2_ball_bound = 2 * metric_r2_bound * metric_r0_bound + 2 * metric_r1_bound**2
        metric_f2_bound = max(
            metric_f2_ball_bound,
            mp.mpf(args.derivative_safety) * metric_sampled_f2_norm,
        )
        metric_data = {
            "metricRowNormSamples": [f(x) for x in metric_norm_samples],
            "metricMidpointFSecondOperatorNorm": f(metric_midpoint_f2_norm),
            "metricSampledFSecondOperatorNorm": f(metric_sampled_f2_norm),
            "metricR3InflatedBound": f(metric_r3_bound),
            "metricR2TaylorBallBound": f(metric_r2_bound),
            "metricR1TaylorBallBound": f(metric_r1_bound),
            "metricR0TaylorBallBound": f(metric_r0_bound),
            "metricFSecondTaylorBallBound": f(metric_f2_ball_bound),
            "metricFSecondSegmentBound": f(metric_f2_bound),
        }

    out = {
        "left": f(left),
        "mid": f(mid),
        "right": f(right),
        "halfWidth": f(half),
        "rowNormSamples": [f(x) for x in norm_samples],
        "midpointFSecondOperatorNorm": f(midpoint_f2_norm),
        "sampledFSecondOperatorNorm": f(sampled_f2_norm),
        "r3InflatedBound": f(r3_bound),
        "r2TaylorBallBound": f(r2_bound),
        "r1TaylorBallBound": f(r1_bound),
        "r0TaylorBallBound": f(r0_bound),
        "fSecondTaylorBallBound": f(f2_ball_bound),
        "fSecondSegmentBound": f(f2_bound),
    }
    if metric_data is not None:
        out.update(metric_data)
    return out


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
    parser.add_argument("--cover-count", type=int, default=6)
    parser.add_argument("--derivative-order", type=int, default=3)
    parser.add_argument("--derivative-safety", default="6")
    parser.add_argument(
        "--frame-json",
        default="trace_frame_finite_gamma_consequence_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_quadrature_interval_consistency_theorem.json",
    )
    parser.add_argument(
        "--reuse-envelope-json",
        default=None,
        help="Reuse a previously computed segment envelope and only recompute theorem bookkeeping.",
    )
    args = parser.parse_args()

    if args.derivative_order < 3:
        raise SystemExit("--derivative-order must be at least 3")

    mp.mp.dps = args.dps
    a = mp.mpf(args.constraint_min)
    b = mp.mpf(args.constraint_max)
    if args.trace_count < 2:
        raise SystemExit("--trace-count must be at least 2")
    h = (b - a) / (args.trace_count - 1)
    trapezoid_factor = (b - a) * h**2 / 12

    frame = load_optional(args.frame_json)
    frame_metric_available = bool(
        frame
        and int(frame.get("basis", -1)) == args.basis
        and int(frame.get("traceCount", -1)) == args.trace_count
        and int(frame.get("activeDimension", -1)) > 0
        and frame.get("centerFrameMatrix")
    )
    metric_inv_sqrt = None
    metric_eigenvalues = None
    if frame_metric_available:
        metric_center = parse_matrix(frame["centerFrameMatrix"])
        metric_inv_sqrt, metric_eigenvalues = inverse_sqrt_spd(metric_center)

    if args.reuse_envelope_json:
        previous = load_optional(args.reuse_envelope_json)
        if not previous:
            raise SystemExit(f"cannot read --reuse-envelope-json {args.reuse_envelope_json}")
        segments = previous["segments"]
        max_f2 = mp.mpf(previous["supFSecondOperatorBound"])
        max_metric_f2 = (
            mp.mpf(previous["supMetricFSecondOperatorBound"])
            if previous.get("supMetricFSecondOperatorBound") is not None
            else None
        )
        active_dimension = int(previous["activeDimension"])
        active_indices = list(previous["activeIndices"])
        active_source_eigs = [mp.mpf(x) for x in previous["activeSourceEigenvalues"]]
        print("Trace quadrature interval consistency theorem")
        print(f"  reused envelope: {args.reuse_envelope_json}")
    else:
        active = build_active_space(args)
        if len(active["activeIndices"]) == 0:
            raise SystemExit("empty active space")

        cover_step = (b - a) / args.cover_count
        segments = []
        max_f2 = mp.mpf("0")
        max_metric_f2 = mp.mpf("0") if metric_inv_sqrt is not None else None
        active_dimension = len(active["activeIndices"])
        active_indices = active["activeIndices"]
        active_source_eigs = active["activeSourceEigenvalues"]
        print("Trace quadrature interval consistency theorem")
        print(f"  basis/traces: {args.basis}/{args.trace_count}")
        print(f"  cover intervals: {args.cover_count}")
        print("  segment       f2_bound")
        for idx in range(args.cover_count):
            left = a + idx * cover_step
            right = a + (idx + 1) * cover_step
            seg = segment_bound(args, active, left, right, metric_inv_sqrt)
            segments.append(seg)
            value = mp.mpf(seg["fSecondSegmentBound"])
            max_f2 = max(max_f2, value)
            if max_metric_f2 is not None:
                max_metric_f2 = max(max_metric_f2, mp.mpf(seg["metricFSecondSegmentBound"]))
            print(f"  {idx:3d} [{fmt(left, 6)}, {fmt(right, 6)}] {fmt(value, 10)}", flush=True)

    error_bound = trapezoid_factor * max_f2
    frame_comparable = bool(
        frame
        and int(frame.get("basis", -1)) == args.basis
        and int(frame.get("traceCount", -1)) == args.trace_count
        and int(frame.get("activeDimension", -1)) == active_dimension
    )
    gamma_finite = mp.mpf(frame["gammaFiniteLower"]) if frame_comparable else None
    lambda_min_center = mp.mpf(frame["lambdaMinCenter"]) if frame_comparable else None
    radius_op = mp.mpf(frame["radiusOperatorBoundUsed"]) if frame_comparable else None
    gamma_after_quadrature = gamma_finite - error_bound if gamma_finite is not None else None
    metric_relative_error = trapezoid_factor * max_metric_f2 if max_metric_f2 is not None else None
    gamma_after_metric_quadrature = (
        (1 - metric_relative_error) * lambda_min_center - radius_op
        if metric_relative_error is not None and lambda_min_center is not None and radius_op is not None
        else None
    )
    quadrature_bound_closed = bool(error_bound >= 0 and max_f2 >= 0)
    gamma_absorbs_error = bool(
        (gamma_after_quadrature is not None and gamma_after_quadrature > 0)
        or (gamma_after_metric_quadrature is not None and gamma_after_metric_quadrature > 0)
    )

    def required_trace_count(target: mp.mpf | None):
        if target is None or target <= 0 or max_f2 <= 0:
            return None
        required_h = mp.sqrt(12 * target / ((b - a) * max_f2))
        return int(mp.ceil(1 + (b - a) / required_h))

    required_for_gamma = required_trace_count(gamma_finite)
    required_for_half_gamma = required_trace_count(gamma_finite / 2) if gamma_finite else None

    data = {
        "theoremName": "trace quadrature interval consistency theorem",
        "proofClass": "interval/ball certificate",
        "basis": args.basis,
        "traceCount": args.trace_count,
        "activeDimension": active_dimension,
        "activeIndices": active_indices,
        "activeSourceEigenvalues": [f(x) for x in active_source_eigs],
        "interval": [f(a), f(b)],
        "meshStep": f(h),
        "coverCount": args.cover_count,
        "derivativeOrder": args.derivative_order,
        "derivativeSafety": f(mp.mpf(args.derivative_safety)),
        "trapezoidFactor": f(trapezoid_factor),
        "supFSecondOperatorBound": f(max_f2),
        "supMetricFSecondOperatorBound": f(max_metric_f2) if max_metric_f2 is not None else None,
        "traceQuadratureErrorBound": f(error_bound),
        "metricRelativeTraceQuadratureErrorBound": f(metric_relative_error)
        if metric_relative_error is not None
        else None,
        "frameJson": args.frame_json if frame else None,
        "frameGammaComparable": frame_comparable,
        "frameMetricEigenvalues": [f(x) for x in metric_eigenvalues] if metric_eigenvalues else None,
        "gammaFiniteLower": f(gamma_finite) if gamma_finite is not None else None,
        "lambdaMinCenter": f(lambda_min_center) if lambda_min_center is not None else None,
        "frameRadiusOperatorBound": f(radius_op) if radius_op is not None else None,
        "gammaAfterTraceQuadrature": f(gamma_after_quadrature)
        if gamma_after_quadrature is not None
        else None,
        "gammaAfterMetricTraceQuadrature": f(gamma_after_metric_quadrature)
        if gamma_after_metric_quadrature is not None
        else None,
        "traceQuadratureClosed": quadrature_bound_closed,
        "traceQuadratureAbsorbedByFiniteGamma": gamma_absorbs_error,
        "requiredTraceCountForErrorBelowGamma": required_for_gamma,
        "requiredTraceCountForErrorBelowHalfGamma": required_for_half_gamma,
        "traceQuadratureConsistencyStatus": {
            "label": "active trace-frame composite trapezoid consistency",
            "closed": quadrature_bound_closed,
            "status": "closed" if quadrature_bound_closed else "open",
            "reason": (
                "The exact active trace derivative rows give a Taylor-ball "
                "enclosure of sup ||(R(a)^*R(a))''|| on the trace interval.  "
                "The composite trapezoid formula gives an explicit operator "
                "error bound for Gamma-Gamma_h."
                if quadrature_bound_closed
                else "The Taylor-ball envelope failed to produce a finite nonnegative error bound."
            ),
        },
        "finiteGammaAbsorptionStatus": {
            "label": "current finite gamma absorbs this quadrature error",
            "closed": gamma_absorbs_error,
            "status": "closed" if gamma_absorbs_error else "open",
            "reason": (
                "The trace quadrature error is controlled either absolutely "
                "below the imported finite interval frame lower bound, or "
                "relatively in the finite frame metric with relative error "
                "below one."
                if gamma_absorbs_error
                else (
                    "The consistency estimate is closed, but this coarse "
                    "trace mesh and conservative derivative envelope consume "
                    "the imported finite gamma.  Use a finer certified trace "
                    "mesh or sharpen the F'' envelope before combining it with "
                    "the finite frame floor."
                )
            ),
        },
        "segments": segments,
        "formula": (
            "||Gamma-Gamma_h|| <= (b-a) h^2/12 sup_a ||d_a^2(R(a)^*R(a))||, "
            "with ||d_a^2(R^*R)|| <= 2||R''||||R||+2||R'||^2.  The metric "
            "version applies the same estimate to R(a) C_h^{-1/2}, where "
            "C_h is the certified finite frame center."
        ),
        "scope": (
            "This certifies trace quadrature consistency on the finite active "
            "Galerkin block.  The remaining separate layer is projection/"
            "Galerkin transport to the continuum active space."
        ),
        "nextProofTarget": (
            "Either certify the finite trace frame on a finer trace mesh or "
            "sharpen the F'' envelope enough that gammaFiniteLower exceeds "
            "the quadrature error, then add the projector transport bound."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print(f"  sup ||F''|| bound: {fmt(max_f2, 12)}")
    if max_metric_f2 is not None:
        print(f"  sup metric ||F''|| bound: {fmt(max_metric_f2, 12)}")
    print(f"  trapezoid factor: {fmt(trapezoid_factor, 12)}")
    print(f"  traceQuadratureErrorBound: {fmt(error_bound, 12)}")
    if metric_relative_error is not None:
        print(f"  metricRelativeErrorBound: {fmt(metric_relative_error, 12)}")
    if gamma_finite is not None:
        print(f"  gammaFiniteLower: {fmt(gamma_finite, 12)}")
        print(f"  gammaAfterTraceQuadrature: {fmt(gamma_after_quadrature, 12)}")
        if gamma_after_metric_quadrature is not None:
            print(f"  gammaAfterMetricTraceQuadrature: {fmt(gamma_after_metric_quadrature, 12)}")
    if gamma_finite is not None:
        print(f"  required trace count for error < gamma: {required_for_gamma}")
        print(f"  required trace count for error < gamma/2: {required_for_half_gamma}")
    print(f"  traceQuadratureClosed: {quadrature_bound_closed}")
    print(f"  absorbed by finite gamma: {gamma_absorbs_error}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
