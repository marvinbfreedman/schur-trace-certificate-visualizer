#!/usr/bin/env python3
r"""Scan the off-base exterior trace determinant landscape.

The one-minor certificate proves finite injectivity, but the continuum proof
needs a stable nonvanishing target.  This script enumerates all off-base
d-by-d sampled trace-response minors on the active block and records the top
determinants.  The columns are normalized so the determinant size measures
angle/independence rather than the scale of the active basis vectors.
"""

from __future__ import annotations

import argparse
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
from trace_active_minor_certificate import normalize_columns, submatrix_rows  # noqa: E402
from trace_to_source_kernel_refinement import active_indices, scale_rows, weights_for_centers  # noqa: E402


def trace_args(args, constraints: int):
    out = child_args(args, args.basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


def all_minors(mat, centers, exclude_center, exclude_radius, zero_tol, top_k):
    dim = mat.cols
    allowed = [
        idx for idx in range(mat.rows)
        if abs(centers[idx] - exclude_center) > exclude_radius
    ]
    rows = []
    sign_counts = {"positive": 0, "negative": 0, "nearZero": 0}
    for combo in itertools.combinations(allowed, dim):
        det = mp.det(submatrix_rows(mat, combo))
        if det > zero_tol:
            sign_counts["positive"] += 1
        elif det < -zero_tol:
            sign_counts["negative"] += 1
        else:
            sign_counts["nearZero"] += 1
        rows.append(
            {
                "rowIndices": list(combo),
                "traceCenters": [f(centers[i]) for i in combo],
                "determinant": f(det),
                "absDeterminant": f(abs(det)),
                "minDistanceToS0": f(min(abs(centers[i] - exclude_center) for i in combo)),
            }
        )
    rows.sort(key=lambda row: mp.mpf(row["absDeterminant"]), reverse=True)
    return rows[:top_k], rows, sign_counts, allowed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--trace-count", type=int, default=9)
    parser.add_argument("--exclude-s0-radius", default="0.001")
    parser.add_argument("--top-k", type=int, default=12)
    parser.add_argument("--zero-tol", default="1e-40")
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
    parser.add_argument("--json-out", default="trace_active_minor_landscape.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    max_local_constraints = max(1, args.basis - args.local_constraint_offset)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    base = local_case(args, args.basis, e_derivs)
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol))
    active_basis = columns(base["sourceVecs"], active_idx)
    centers, R_global = trace_matrix(base["polys"], trace_args(args, args.trace_count))
    trace_active = R_global * base["scaledModes"] * active_basis
    weights = weights_for_centers(centers)
    weighted_trace = scale_rows(trace_active, [mp.sqrt(w) for w in weights])
    normalized, column_norms = normalize_columns(weighted_trace)
    top, all_rows, sign_counts, allowed = all_minors(
        normalized,
        centers,
        mp.mpf(args.s0),
        mp.mpf(args.exclude_s0_radius),
        mp.mpf(args.zero_tol),
        args.top_k,
    )

    frame_vals = mp.eigsy((weighted_trace.T * weighted_trace + weighted_trace.T * weighted_trace) / 2, eigvals_only=True)
    best = top[0] if top else None
    kth = top[-1] if top else None
    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basis": args.basis,
        "traceCount": args.trace_count,
        "excludeS0Radius": f(mp.mpf(args.exclude_s0_radius)),
        "allowedTraceCenters": [f(centers[i]) for i in allowed],
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "weightedColumnNorms": [f(x) for x in column_norms],
        "weightedFrameEigenvalues": [f(x) for x in frame_vals],
        "minorCount": len(all_rows),
        "signCounts": sign_counts,
        "bestAbsDeterminant": best["absDeterminant"] if best else None,
        "topKthAbsDeterminant": kth["absDeterminant"] if kth else None,
        "topMinors": top,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "interpretation": (
            "All listed determinants use trace points outside the excluded "
            "s0-neighborhood.  Persistent nonzero top minors are the finite "
            "shadow of the exterior trace nonvanishing theorem."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print(
        f"Trace active minor landscape basis={args.basis} trace_count={args.trace_count} "
        f"active_dim={len(active_idx)} minors={len(all_rows)}"
    )
    print(f"  frame_min={fmt(frame_vals[0], 10) if len(frame_vals) else '0'}")
    print(f"  sign counts={sign_counts}")
    print("  rank abs_det centers")
    for i, row in enumerate(top, start=1):
        print(
            f"  {i:4d} {fmt(mp.mpf(row['absDeterminant']), 10):>12} "
            f"{row['traceCenters']}",
            flush=True,
        )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
