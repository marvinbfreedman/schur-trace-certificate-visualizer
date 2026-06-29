#!/usr/bin/env python3
r"""Trace resolution test for the Mellin boundary prefix.

The exact atom identity is

    X_i(z)=B_i(s,z)+T_i(s,z),

where T_i is the diagonal Volterra tail and B_i is an incomplete-gamma prefix
from the fixed lower endpoint to the moving base point s.  The diagonal tail
already lives in the Volterra atom system.  This script tests the next needed
statement:

    boundary-prefix source rows are endpoint-trace rows.

Finite model.  For each de Branges sample z and branch shift z -> z +/- i omega,
build the Legendre coefficient row

    E_B(z)_k = int_0^L B(s,z) p_k(s) ds,

then solve

    E_B(z) ~= c_z^T R,

where R is the sampled endpoint trace matrix

    (R f)_j = Lambda_{a_j}(f).

If this row-span residual is small in the correct quotient/energy norm, then
the boundary prefix vanishes on ker R and the remaining transport is the
diagonal Volterra tail from the exact Mellin identity.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import default_z_nodes  # noqa: E402
from klm_debranges_trace_map_constructor import ctranspose, to_float  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    gram_matrix,
    legendre_quadrature,
    poly_value,
    shifted_legendre_polys,
    trace_matrix,
)
from xi_mellin_convolution_boundary_identity import atom_prefix  # noqa: E402
from xi_mellin_volterra_mode_match import atom_records  # noqa: E402


def local_trace_args(args: argparse.Namespace, constraints: int | None = None):
    return SimpleNamespace(
        constraints=args.constraints if constraints is None else constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
    )


def branch_z_values(args):
    omega = mp.mpf(args.omega)
    out = []
    for idx, z in enumerate(default_z_nodes()):
        out.append((idx, "plus", mp.mpc(z) + mp.j * omega))
        out.append((idx, "minus", mp.mpc(z) - mp.j * omega))
    return out


def prefix_sum(records, s: mp.mpf, z, fourier_scale: mp.mpf, amp_kind: str) -> mp.mpc:
    total = mp.fsum(atom_prefix(record, s, z, fourier_scale) for record in records)
    if amp_kind == "positive":
        return total
    if amp_kind == "even":
        return total + mp.fsum(atom_prefix(record, s, -mp.mpc(z), fourier_scale) for record in records)
    raise ValueError(amp_kind)


def boundary_prefix_row(args, polys, z) -> mp.matrix:
    records = atom_records(args.kind)
    L = mp.mpf(args.L)
    s_pts, s_wts = legendre_quadrature(L, args.quad)
    row = mp.matrix(1, len(polys))
    scale = mp.mpf(args.fourier_scale)
    norm = mp.sqrt(2 * mp.pi) if args.hardy_normalization else mp.mpf("1")
    for s, w in zip(s_pts, s_wts):
        val = prefix_sum(records, s, z, scale, args.amp_kind) / norm
        for k, poly in enumerate(polys):
            row[0, k] += w * val * poly_value(poly, s)
    return row


def solve_row_span(R: mp.matrix, row: mp.matrix, reg: mp.mpf):
    # Minimize ||R^T c - row^T||_2.
    gram = R * R.T
    rhs = R * row.T
    diag_max = max([abs(gram[i, i]) for i in range(gram.rows)] + [mp.mpf("1")])
    lhs = gram.copy()
    for i in range(lhs.rows):
        lhs[i, i] += reg * diag_max
    coeff = mp.lu_solve(lhs, rhs) if lhs.rows else mp.matrix(0, 1)
    recon = coeff.T * R if lhs.rows else mp.matrix(1, row.cols)
    resid = row - recon
    return coeff, recon, resid


def split_trace_nullspace(R: mp.matrix, rank_tol_text: str):
    gram = R.T * R
    vals, vecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    vmax = max(abs(v) for v in vals) if len(vals) else mp.mpf("0")
    tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), vmax)
    n_idx = [i for i, val in enumerate(vals) if val <= tol]
    u_idx = [i for i, val in enumerate(vals) if val > tol]
    return columns(vecs, n_idx), columns(vecs, u_idx), vals, tol


def weighted_null_residual(row: mp.matrix, N: mp.matrix, K: mp.matrix, psd_tol: mp.mpf):
    if N.cols == 0:
        return mp.mpf("0"), mp.mpf("0"), 0
    A = N.T * K * N
    vals, vecs = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > psd_tol]
    if not keep:
        return mp.mpf("0"), mp.mpf("0"), 0
    V = columns(vecs, keep)
    max_ratio = mp.mpf("0")
    sum_ratio = mp.mpf("0")
    for local_idx, eig_idx in enumerate(keep):
        mode = N * V[:, local_idx]
        amp = (row * mode)[0]
        ratio = abs(amp) / mp.sqrt(vals[eig_idx])
        max_ratio = max(max_ratio, ratio)
        sum_ratio += ratio * ratio
    return max_ratio, mp.sqrt(sum_ratio), len(keep)


def compute(args):
    L = mp.mpf(args.L)
    polys = shifted_legendre_polys(args.basis, L)
    _centers, R = trace_matrix(polys, local_trace_args(args))
    K, _kpolys = gram_matrix(args, mp.mpf(args.omega), L)
    N, U, rvals, rank_tol = split_trace_nullspace(R, args.rank_tol)
    rows = []
    max_rel = mp.mpf("0")
    max_abs = mp.mpf("0")
    max_null = mp.mpf("0")
    max_null_l2 = mp.mpf("0")
    for node_index, branch, z in branch_z_values(args):
        row = boundary_prefix_row(args, polys, z)
        coeff, recon, resid = solve_row_span(R, row, mp.mpf(args.trace_reg))
        row_norm = mp.sqrt(mp.fsum(abs(row[0, i]) ** 2 for i in range(row.cols)))
        resid_norm = mp.sqrt(mp.fsum(abs(resid[0, i]) ** 2 for i in range(resid.cols)))
        rel = resid_norm / max(row_norm, mp.mpf("1e-300"))
        null_max, null_l2, positive_modes = weighted_null_residual(
            row, N, K, mp.mpf(args.psd_tol)
        )
        max_rel = max(max_rel, rel)
        max_abs = max(max_abs, resid_norm)
        max_null = max(max_null, null_max)
        max_null_l2 = max(max_null_l2, null_l2)
        rows.append(
            {
                "nodeIndex": node_index,
                "branch": branch,
                "z": {"real": to_float(mp.re(z)), "imag": to_float(mp.im(z))},
                "rowNorm": to_float(row_norm),
                "rowSpanResidualNorm": to_float(resid_norm),
                "rowSpanResidualRelative": to_float(rel),
                "traceCoefficientNorm": to_float(mp.sqrt(mp.fsum(abs(coeff[i]) ** 2 for i in range(coeff.rows)))),
                "nullEnergyOpMax": to_float(null_max),
                "nullEnergyOpL2": to_float(null_l2),
                "positiveNullModes": positive_modes,
            }
        )
    return {
        "traceRank": U.cols,
        "traceNullity": N.cols,
        "traceRankTolerance": to_float(rank_tol),
        "maxRowSpanResidualRelative": to_float(max_rel),
        "maxRowSpanResidualNorm": to_float(max_abs),
        "maxNullEnergyOpMax": to_float(max_null),
        "maxNullEnergyOpL2": to_float(max_null_l2),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="full")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=10)
    parser.add_argument("--quad", type=int, default=18)
    parser.add_argument("--constraints", type=int, default=7)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=20)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-24")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--trace-reg", default="1e-40")
    parser.add_argument("--fourier-scale", default="0.5")
    parser.add_argument("--amp-kind", choices=["positive", "even"], default="positive")
    parser.add_argument("--hardy-normalization", action="store_true", default=True)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="xi_boundary_prefix_trace_resolution.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    result = compute(args)
    data = {
        "theoremName": "Mellin boundary-prefix endpoint trace resolution",
        "kind": args.kind,
        "omega": float(args.omega),
        "L": float(args.L),
        "basis": args.basis,
        "constraints": args.constraints,
        "ampKind": args.amp_kind,
        "fourierScale": float(args.fourier_scale),
        **result,
        "boundaryPrefixTraceResolved": result["maxRowSpanResidualRelative"] < 1e-8,
        "boundaryPrefixKilledOnTraceNullspace": result["maxNullEnergyOpMax"] < 1e-8,
        "diagnosis": (
            "This tests whether the incomplete-gamma boundary prefix lies in the "
            "finite sampled endpoint-trace row span.  The null-energy residual is "
            "the relevant quotient diagnostic: if it vanishes, the prefix is "
            "zero on ker R even if the Euclidean row-span residual is not tiny."
        ),
        "nextProofTarget": (
            "If the boundary prefix is not killed on the trace nullspace, enlarge "
            "the trace family to include the Mellin boundary primitive or derive "
            "the missing endpoint concomitant that pairs with B_i(s,z)."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Mellin boundary-prefix endpoint trace resolution")
    print(
        f"  kind={args.kind} omega={args.omega} basis={args.basis} "
        f"constraints={args.constraints} nullity={result['traceNullity']}"
    )
    print(f"  max row-span residual relative: {result['maxRowSpanResidualRelative']:.3e}")
    print(f"  max null-energy op max: {result['maxNullEnergyOpMax']:.3e}")
    print(f"  max null-energy op l2: {result['maxNullEnergyOpL2']:.3e}")
    print(f"  trace resolved: {data['boundaryPrefixTraceResolved']}")
    print(f"  killed on trace nullspace: {data['boundaryPrefixKilledOnTraceNullspace']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
