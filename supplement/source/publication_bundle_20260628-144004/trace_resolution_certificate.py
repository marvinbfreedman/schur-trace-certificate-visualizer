#!/usr/bin/env python3
r"""Finite certificate for the endpoint trace-resolution lemma.

For F(a)=Lambda_a(f), sampled on nodes a_j, the elementary mesh estimate is

    sup_a |F(a)| <= max_j |F(a_j)| + delta sup_a |F'(a)|,

where delta is the mesh fill distance.  In a spectral block E_M of
A=K|ker(R), this becomes the operator estimate

    ||T_dense||_{A->l_inf}
      <= ||T_sample||_{A->l_inf} + delta ||dT/da||_{A->l_inf}.

This script computes those three constants for low/mid A-spectral blocks and
also reports jet constants

    sup_a |f^(k)(a)| <= J_{k,M} ||f||_A,

which are the finite-dimensional energy-to-jet estimates needed to bound
d/da Lambda_a(f) analytically.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import (  # noqa: E402
    columns,
    gram_matrix,
    poly_derivative_value,
    trace_matrix,
)


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def f(x):
    return float(x)


def fmt(x, digits=8):
    return mp.nstr(x, digits)


def make_args(base, constraints):
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=base.basis,
        quad=base.quad,
        laguerre=base.laguerre,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
        dps=base.dps,
    )


def split_spaces(R, rank_tol_text):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx), len(u_idx), len(n_idx)


def centers_for(args):
    lo = mp.mpf(args.constraint_min)
    hi = mp.mpf(args.constraint_max)
    if args.constraints == 1:
        return [(lo + hi) / 2]
    step = (hi - lo) / (args.constraints - 1)
    return [lo + i * step for i in range(args.constraints)]


def fill_distance(dense_centers, sample_centers):
    return max(min(abs(a - b) for b in sample_centers) for a in dense_centers)


def weighted_modes(modes, lambdas, cutoff):
    out = mp.matrix(modes.rows, cutoff)
    for j in range(cutoff):
        scale = 1 / mp.sqrt(lambdas[j])
        for i in range(modes.rows):
            out[i, j] = modes[i, j] * scale
    return out


def row_op_norm(row, wmodes):
    vals = []
    for j in range(wmodes.cols):
        vals.append(mp.fsum(row[0, i] * wmodes[i, j] for i in range(wmodes.rows)))
    return mp.sqrt(mp.fsum(v * v for v in vals))


def max_trace_op_norm(R, wmodes):
    if R.rows == 0 or wmodes.cols == 0:
        return mp.mpf("0")
    vals = []
    for i in range(R.rows):
        row = mp.matrix(1, R.cols)
        for j in range(R.cols):
            row[0, j] = R[i, j]
        vals.append(row_op_norm(row, wmodes))
    return max(vals)


def derivative_rows(R, centers):
    out = mp.matrix(R.rows, R.cols)
    for i in range(R.rows):
        if i == 0:
            h = centers[1] - centers[0]
            for j in range(R.cols):
                out[i, j] = (R[1, j] - R[0, j]) / h
        elif i == R.rows - 1:
            h = centers[-1] - centers[-2]
            for j in range(R.cols):
                out[i, j] = (R[-1, j] - R[-2, j]) / h
        else:
            h = centers[i + 1] - centers[i - 1]
            for j in range(R.cols):
                out[i, j] = (R[i + 1, j] - R[i - 1, j]) / h
    return out


def jet_row(polys, x, deriv):
    row = mp.matrix(1, len(polys))
    for i, poly in enumerate(polys):
        row[0, i] = poly_derivative_value(poly, x, deriv)
    return row


def jet_constants(polys, dense_centers, wmodes, max_deriv):
    constants = []
    for deriv in range(max_deriv + 1):
        best = mp.mpf("0")
        for x in dense_centers:
            best = max(best, row_op_norm(jet_row(polys, x, deriv), wmodes))
        constants.append(best)
    return constants


def compute_case(base, K, polys, dense_R, dense_centers, constraints, cutoff):
    args = make_args(base, constraints)
    sample_centers, R = trace_matrix(polys, args)
    N, U, rank, nullity = split_spaces(R, args.rank_tol)
    if N.cols == 0:
        return None
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    if not keep:
        return None
    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    modes = N * avecs
    cutoff = min(cutoff, len(avals))
    if cutoff == 0:
        return None
    wmodes = weighted_modes(modes, avals, cutoff)
    dR = derivative_rows(dense_R, dense_centers)

    sample_op = max_trace_op_norm(R, wmodes)
    dense_op = max_trace_op_norm(dense_R, wmodes)
    deriv_op = max_trace_op_norm(dR, wmodes)
    delta = fill_distance(dense_centers, sample_centers)
    bound = sample_op + delta * deriv_op
    jets = jet_constants(polys, dense_centers, wmodes, base.max_jet_deriv)
    return {
        "constraints": constraints,
        "rank": rank,
        "nullity": nullity,
        "cutoff": cutoff,
        "positiveModes": len(avals),
        "fillDistance": f(delta),
        "sampleTraceOp": f(sample_op),
        "denseTraceOp": f(dense_op),
        "traceDerivativeOp": f(deriv_op),
        "meshBound": f(bound),
        "denseOverBound": f(dense_op / bound) if bound else 0.0,
        "jetConstants": [f(x) for x in jets],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=16)
    parser.add_argument("--quad", type=int, default=22)
    parser.add_argument("--constraints", default="8 10 12")
    parser.add_argument("--cutoff", type=int, default=7)
    parser.add_argument("--dense-constraints", type=int, default=33)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--max-jet-deriv", type=int, default=10)
    parser.add_argument("--json-out", default="trace_resolution_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    length = mp.mpf(args.L)
    K, polys = gram_matrix(args, mp.mpf(args.omega), length)
    dense_args = make_args(args, args.dense_constraints)
    dense_centers, dense_R = trace_matrix(polys, dense_args)

    rows = []
    print(
        f"Trace resolution certificate model={args.model} kind={args.kind} "
        f"basis={args.basis} cutoff={args.cutoff}",
        flush=True,
    )
    print(
        "  cons cutoff delta      sample_op  dense_op   dtrace_op  bound      ratio",
        flush=True,
    )
    for constraints in parse_ints(args.constraints):
        row = compute_case(args, K, polys, dense_R, dense_centers, constraints, args.cutoff)
        if row is None:
            continue
        rows.append(row)
        print(
            f"  {constraints:4d} {row['cutoff']:6d} {fmt(row['fillDistance']):>10} "
            f"{fmt(row['sampleTraceOp']):>10} {fmt(row['denseTraceOp']):>10} "
            f"{fmt(row['traceDerivativeOp']):>10} {fmt(row['meshBound']):>10} "
            f"{fmt(row['denseOverBound']):>8}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(length),
        "basis": args.basis,
        "denseConstraints": args.dense_constraints,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
