#!/usr/bin/env python3
r"""Certificate for finite active trace injectivity by trace-response minors.

On the source-active block H_delta, the moving trace map is represented by a
small matrix whose columns are the trace-response functions

    a -> Lambda_a(v_k).

The weighted frame bound proves injectivity through T^*T.  This script records
the sharper finite-dimensional fact: a single nonzero d-by-d minor of the
sampled trace-response matrix already proves injectivity on the d-dimensional
active block.  In the continuum theorem this is the Chebyshev/unique-continuity
object to prove.
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
from global_trace_observability_gap import child_args, f, fmt, parse_ints  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns, trace_matrix  # noqa: E402
from trace_to_source_kernel_refinement import active_indices, scale_rows, weights_for_centers  # noqa: E402


def trace_args(args, basis: int, constraints: int):
    out = child_args(args, basis, constraints)
    out.constraints = constraints
    return make_qargs(out)


def submatrix_rows(mat, row_indices):
    out = mp.matrix(len(row_indices), mat.cols)
    for i, row in enumerate(row_indices):
        for j in range(mat.cols):
            out[i, j] = mat[row, j]
    return out


def column_norms(mat):
    return [
        mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows)))
        for j in range(mat.cols)
    ]


def normalize_columns(mat):
    norms = column_norms(mat)
    out = mp.matrix(mat.rows, mat.cols)
    for j, norm in enumerate(norms):
        scale = 1 / norm if norm else mp.mpf("0")
        for i in range(mat.rows):
            out[i, j] = mat[i, j] * scale
    return out, norms


def find_best_minor(mat, centers, exclude_center=None, exclude_radius=mp.mpf("0")):
    dim = mat.cols
    if dim == 0:
        return None
    allowed = list(range(mat.rows))
    if exclude_center is not None and exclude_radius > 0:
        allowed = [
            idx for idx in allowed
            if abs(centers[idx] - exclude_center) > exclude_radius
        ]
    if len(allowed) < dim:
        return None
    best = None
    for rows in itertools.combinations(allowed, dim):
        sub = submatrix_rows(mat, rows)
        det = mp.det(sub)
        item = {
            "rowIndices": list(rows),
            "traceCenters": [f(centers[i]) for i in rows],
            "determinant": f(det),
            "absDeterminant": f(abs(det)),
        }
        if best is None or abs(det) > mp.mpf(best["absDeterminant"]):
            best = item
    return best


def frame_eigenvalues(mat):
    if mat.cols == 0:
        return []
    vals = mp.eigsy((mat.T * mat + mat.T * mat) / 2, eigvals_only=True)
    return [f(v) for v in vals]


def run_case(args, basis: int, trace_count: int, e_derivs):
    base = local_case(args, basis, e_derivs)
    active_idx, active_floor = active_indices(base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol))
    active_basis = columns(base["sourceVecs"], active_idx)
    qargs = trace_args(args, basis, trace_count)
    centers, R_global = trace_matrix(base["polys"], qargs)
    trace_active = R_global * base["scaledModes"] * active_basis
    weights = weights_for_centers(centers)
    weighted_trace = scale_rows(trace_active, [mp.sqrt(w) for w in weights])
    normalized, norms = normalize_columns(weighted_trace)
    best = find_best_minor(
        normalized,
        centers,
        mp.mpf(args.s0),
        mp.mpf(args.exclude_s0_radius),
    )
    frame_vals = frame_eigenvalues(weighted_trace)
    normalized_frame_vals = frame_eigenvalues(normalized)
    return {
        "basis": basis,
        "traceCount": trace_count,
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "activeDim": len(active_idx),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "traceCenters": [f(x) for x in centers],
        "weightedColumnNorms": [f(x) for x in norms],
        "weightedFrameEigenvalues": frame_vals,
        "normalizedFrameEigenvalues": normalized_frame_vals,
        "bestNormalizedMinor": best,
        "interpretation": (
            "A nonzero bestNormalizedMinor proves finite sampled trace "
            "injectivity on this active block.  The normalized minor removes "
            "irrelevant scaling of the active basis vectors."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="22")
    parser.add_argument("--trace-count", type=int, default=9)
    parser.add_argument("--exclude-s0-radius", default="0")
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
    parser.add_argument("--json-out", default="trace_active_minor_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(1, max(bases) - args.local_constraint_offset)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)

    print(f"Trace active minor certificate bases={bases} trace_count={args.trace_count}")
    print("  basis active_dim frame_min best_abs_minor centers")
    rows = []
    for basis in bases:
        row = run_case(args, basis, args.trace_count, e_derivs)
        rows.append(row)
        frame_min = mp.mpf(row["weightedFrameEigenvalues"][0]) if row["weightedFrameEigenvalues"] else mp.mpf("0")
        minor = row["bestNormalizedMinor"] or {}
        print(
            f"  {basis:5d} {row['activeDim']:10d} {fmt(frame_min, 8):>10} "
            f"{fmt(mp.mpf(minor.get('absDeterminant', 0)), 8):>14} "
            f"{minor.get('traceCenters', [])}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "basisList": bases,
        "traceCount": args.trace_count,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8AtS0": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
