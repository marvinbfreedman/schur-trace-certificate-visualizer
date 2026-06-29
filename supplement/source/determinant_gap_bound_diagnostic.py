#!/usr/bin/env python3
r"""Finite diagnostic for the determinant-gap inequality.

The abstract continuum step is

    ||M_N - M|| <= C alpha_N < sigma_min(M_N),

where C is the A-dual norm of the chosen derivative rows and alpha_N is the
active spectral-subspace graph error.  The derivative-rank scans also report a
column-normalized singular value; that number is useful for angles, but it is
not in the same units as the raw inequality above.

This script computes the finite raw quantities for a fixed off-base tuple:

    a0=.545, q=(7,8)

on the stable two-dimensional source-active band, plus a Galerkin subspace
drift diagnostic in a common parent polynomial space.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import (  # noqa: E402
    exact_trace_derivatives,
    tower_row_on_polys,
)
from quotient_factorization_mp import columns, gram_matrix  # noqa: E402
from trace_active_derivative_rank import derivative_matrix, normalized_columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def parse_orders(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def eigvals_gram(mat):
    return mp.eigsy((mat.T * mat + mat.T * mat) / 2, eigvals_only=True)


def submatrix_rows(mat, rows):
    out = mp.matrix(len(rows), mat.cols)
    for i, row in enumerate(rows):
        for j in range(mat.cols):
            out[i, j] = mat[row, j]
    return out


def pad_rows(mat, rows: int):
    out = mp.matrix(rows, mat.cols)
    for i in range(mat.rows):
        for j in range(mat.cols):
            out[i, j] = mat[i, j]
    return out


def sym(mat):
    return (mat + mat.T) / 2


def k_orthonormalize(mat, K, tol):
    gram = sym(mat.T * K * mat)
    vals, vecs = mp.eigsy(gram, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > tol]
    if not keep:
        return mp.matrix(mat.rows, 0), []
    vkeep = columns(vecs, keep)
    scales = mp.matrix(len(keep))
    for j, idx in enumerate(keep):
        scales[j, j] = 1 / mp.sqrt(vals[idx])
    return mat * vkeep * scales, [vals[i] for i in keep]


def principal_angle_drift(left, right, K, tol):
    ql, lvals = k_orthonormalize(left, K, tol)
    qr, rvals = k_orthonormalize(right, K, tol)
    if ql.cols == 0 or qr.cols == 0:
        return None
    cross = ql.T * K * qr
    vals = eigvals_gram(cross)
    svals = [mp.sqrt(max(mp.mpf("0"), val)) for val in vals]
    smin = min(svals) if svals else mp.mpf("0")
    smax = max(svals) if svals else mp.mpf("0")
    projection_gap = mp.sqrt(max(mp.mpf("0"), 1 - smin**2))
    polar_alpha = mp.sqrt(max(mp.mpf("0"), 2 - 2 * smin))
    return {
        "leftDim": ql.cols,
        "rightDim": qr.cols,
        "leftGramEigenvalues": [f(x) for x in lvals],
        "rightGramEigenvalues": [f(x) for x in rvals],
        "singularValues": [f(x) for x in svals],
        "minCosine": f(smin),
        "maxCosine": f(smax),
        "projectionGapEstimate": f(projection_gap),
        "polarAlphaEstimate": f(polar_alpha),
    }


def basis_case(args, basis: int, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(
        base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol)
    )
    active_basis = columns(base["sourceVecs"], active_idx)
    active_modes = base["scaledModes"] * active_basis

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, e_derivs_point, _lam_derivs = exact_trace_derivatives(
        point_args, max(args.orders)
    )
    full_matrix = derivative_matrix(
        point_args, base, active_modes, mp.mpf(args.point), max(args.orders)
    )
    raw = submatrix_rows(full_matrix, args.orders)
    normalized, column_norms = normalized_columns(raw)
    raw_vals = eigvals_gram(raw)
    norm_vals = eigvals_gram(normalized)
    raw_smin = mp.sqrt(max(mp.mpf("0"), raw_vals[0])) if raw_vals else mp.mpf("0")
    normalized_smin = mp.sqrt(max(mp.mpf("0"), norm_vals[0])) if norm_vals else mp.mpf("0")

    row_norm2_sum = mp.mpf("0")
    row_norms = []
    for q in args.orders:
        row = tower_row_on_polys(
            base["polys"], mp.mpf(args.point), e_derivs_point, q, args.jet_order
        )
        on_high = row * base["scaledModes"]
        norm = mp.sqrt(mp.fsum(abs(on_high[0, j]) ** 2 for j in range(on_high.cols)))
        row_norms.append(norm)
        row_norm2_sum += norm**2
    C_high = mp.sqrt(row_norm2_sum)
    allowed_alpha_raw = raw_smin / C_high if C_high else mp.mpf("0")
    return {
        "basis": basis,
        "activeDim": len(active_idx),
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "activeSourceEigenvalues": [f(base["sourceVals"][i]) for i in active_idx],
        "rawSingularValuesSquared": [f(x) for x in raw_vals],
        "normalizedSingularValuesSquared": [f(x) for x in norm_vals],
        "rawMinSingularValue": f(raw_smin),
        "normalizedMinSingularValue": f(normalized_smin),
        "columnNorms": [f(x) for x in column_norms],
        "rowNormsOnHighBlock": [f(x) for x in row_norms],
        "COnHighBlock": f(C_high),
        "allowedAlphaRaw": f(allowed_alpha_raw),
        "activeModes": active_modes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20 22 24")
    parser.add_argument("--basis", type=int, default=24)
    parser.add_argument("--max-deriv", type=int, default=8)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-6")
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
    parser.add_argument("--json-out", default="determinant_gap_bound_diagnostic.json")
    raw_args = parser.parse_args()

    raw_args.orders = parse_orders(raw_args.orders)
    mp.mp.dps = raw_args.dps
    bases = parse_ints(raw_args.bases)
    max_basis = max(bases)
    max_local_constraints = max(max(1, b - raw_args.local_constraint_offset) for b in bases)
    max_q = max(raw_args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(raw_args, max_q)

    print(
        f"Determinant gap bound diagnostic point={raw_args.point} "
        f"orders={raw_args.orders} bases={bases}",
        flush=True,
    )
    rows = []
    active_modes_by_basis = {}
    for basis in bases:
        row = basis_case(raw_args, basis, e_derivs_at_s0)
        active_modes_by_basis[basis] = row.pop("activeModes")
        rows.append(row)
        print(
            "  basis {basis:2d}, dim {dim}: raw_smin={raw}, C={C}, "
            "allowed_alpha={alpha}, norm_smin={norm}".format(
                basis=basis,
                dim=row["activeDim"],
                raw=fmt(mp.mpf(row["rawMinSingularValue"]), 8),
                C=fmt(mp.mpf(row["COnHighBlock"]), 8),
                alpha=fmt(mp.mpf(row["allowedAlphaRaw"]), 8),
                norm=fmt(mp.mpf(row["normalizedMinSingularValue"]), 8),
            ),
            flush=True,
        )

    parent_args = SimpleNamespace(**vars(raw_args))
    parent_args.basis = max_basis
    parent_args.constraints = max(1, max_basis - raw_args.local_constraint_offset)
    K_parent, _polys = gram_matrix(
        make_qargs(parent_args), mp.mpf(raw_args.omega), mp.mpf(raw_args.L)
    )
    drifts = []
    for left, right in zip(bases[:-1], bases[1:]):
        lpad = pad_rows(active_modes_by_basis[left], max_basis)
        rpad = pad_rows(active_modes_by_basis[right], max_basis)
        drift = principal_angle_drift(lpad, rpad, K_parent, mp.mpf(raw_args.psd_tol))
        if drift is not None:
            drift["leftBasis"] = left
            drift["rightBasis"] = right
            drifts.append(drift)
            print(
                f"  drift {left}->{right}: min_cos={fmt(mp.mpf(drift['minCosine']), 8)} "
                f"polar_alpha={fmt(mp.mpf(drift['polarAlphaEstimate']), 8)}",
                flush=True,
            )

    min_normalized = min(mp.mpf(row["normalizedMinSingularValue"]) for row in rows)
    min_allowed_raw = min(mp.mpf(row["allowedAlphaRaw"]) for row in rows)
    max_drift = max(
        [mp.mpf(row["polarAlphaEstimate"]) for row in drifts] + [mp.mpf("0")]
    )
    data = {
        "omega": f(mp.mpf(raw_args.omega)),
        "L": f(mp.mpf(raw_args.L)),
        "s0": f(mp.mpf(raw_args.s0)),
        "point": f(mp.mpf(raw_args.point)),
        "orders": raw_args.orders,
        "bases": bases,
        "activeTol": f(mp.mpf(raw_args.active_tol)),
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "minNormalizedSingularValue": f(min_normalized),
        "minAllowedAlphaRaw": f(min_allowed_raw),
        "maxConsecutivePolarAlphaEstimate": f(max_drift),
        "rows": rows,
        "drifts": drifts,
        "interpretation": (
            "The normalized singular margin is not the same quantity as the "
            "raw perturbation bound C alpha.  The raw finite theorem would "
            "need alpha below allowedAlphaRaw, or a reformulated normalized "
            "perturbation estimate with explicit column-scale control."
        ),
    }
    Path(raw_args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {raw_args.json_out}")


if __name__ == "__main__":
    main()
