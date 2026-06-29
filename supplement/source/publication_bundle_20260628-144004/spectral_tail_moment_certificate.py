#!/usr/bin/env python3
r"""Moment certificate for spectral high-tail decay.

For A=K|ker(R) and source block B, the spectral theorem gives for q>0:

    C_tail(Lambda)
      = || A^{-1/2} 1_{A>=Lambda} B ||^2
      <= Lambda^{-q} || A^{(q-1)/2} B ||^2.

This script computes the source moments

    S_q = || A^{(q-1)/2} B ||^2

and compares the resulting bounds with the observed spectral tails.  This is
the finite certificate for the high-frequency Hardy tail theorem.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from high_frequency_tail_refinement import constraints_for  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, max_eig_or_zero, trace_matrix  # noqa: E402


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def parse_mpf_list(text):
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def f(x):
    return float(x)


def fmt(x, digits=8):
    return mp.nstr(x, digits)


def make_args(base, basis, constraints):
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=base.quad_base + base.quad_step * basis,
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


def row_outer(row, scale):
    out = mp.matrix(row.cols)
    for i in range(row.cols):
        for j in range(row.cols):
            out[i, j] = scale * row[0, i] * row[0, j]
    return out


def weighted_source_constant(beig, avals, power, threshold=None):
    if beig.cols == 0:
        return mp.mpf("0")
    mat = mp.matrix(beig.cols)
    for j, lam in enumerate(avals):
        if threshold is not None and lam < threshold:
            continue
        row = mp.matrix(1, beig.cols)
        for col in range(beig.cols):
            row[0, col] = beig[j, col]
        mat += row_outer(row, lam ** power)
    return max_eig_or_zero(mat)


def cutoff_for_threshold(avals, threshold):
    for idx, lam in enumerate(avals):
        if lam >= threshold:
            return idx
    return len(avals)


def compute_case(base, basis):
    constraints = constraints_for(base, basis)
    args = make_args(base, basis, constraints)
    K, polys = gram_matrix(args, mp.mpf(args.omega), mp.mpf(args.L))
    _, R = trace_matrix(polys, args)
    N, U, rank, nullity = split_spaces(R, args.rank_tol)
    if N.cols == 0 or U.cols == 0:
        return None

    A = N.T * K * N
    B = N.T * K * U
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    if not keep:
        return None
    avals = [avals_all[i] for i in keep]
    avecs = columns(avecs_all, keep)
    beig = avecs.T * B

    full = weighted_source_constant(beig, avals, mp.mpf("-1"))
    moments = []
    qs = parse_mpf_list(base.q_values)
    for q in qs:
        sq = weighted_source_constant(beig, avals, q - 1)
        moments.append({"q": q, "S": sq, "relative": sq / full if full else mp.mpf("0")})

    thresholds = []
    for threshold in parse_mpf_list(base.thresholds):
        tail = weighted_source_constant(beig, avals, mp.mpf("-1"), threshold=threshold)
        bounds = []
        for item in moments:
            bound = threshold ** (-item["q"]) * item["S"]
            bounds.append(
                {
                    "q": item["q"],
                    "bound": bound,
                    "boundFrac": bound / full if full else mp.mpf("0"),
                    "efficiency": tail / bound if bound else mp.mpf("0"),
                }
            )
        thresholds.append(
            {
                "threshold": threshold,
                "cutoff": cutoff_for_threshold(avals, threshold),
                "tail": tail,
                "tailFrac": tail / full if full else mp.mpf("0"),
                "bounds": bounds,
            }
        )

    return {
        "basis": basis,
        "constraints": constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "gamma2": full,
        "lambdaMin": avals[0],
        "lambdaMax": avals[-1],
        "moments": moments,
        "thresholds": thresholds,
    }


def serial(row):
    return {
        "basis": row["basis"],
        "constraints": row["constraints"],
        "rank": row["rank"],
        "nullity": row["nullity"],
        "positiveModes": row["positiveModes"],
        "gamma2": f(row["gamma2"]),
        "lambdaMin": f(row["lambdaMin"]),
        "lambdaMax": f(row["lambdaMax"]),
        "moments": [
            {
                "q": f(item["q"]),
                "S": f(item["S"]),
                "relative": f(item["relative"]),
            }
            for item in row["moments"]
        ],
        "thresholds": [
            {
                "threshold": f(item["threshold"]),
                "cutoff": item["cutoff"],
                "tail": f(item["tail"]),
                "tailFrac": f(item["tailFrac"]),
                "bounds": [
                    {
                        "q": f(bound["q"]),
                        "bound": f(bound["bound"]),
                        "boundFrac": f(bound["boundFrac"]),
                        "efficiency": f(bound["efficiency"]),
                    }
                    for bound in item["bounds"]
                ],
            }
            for item in row["thresholds"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="12 14 16 18")
    parser.add_argument("--constraint-rule", choices=["ratio", "target-nullity", "offset"], default="ratio")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--target-nullity", type=int, default=6)
    parser.add_argument("--constraint-offset", type=int, default=6)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--quad-base", type=int, default=6)
    parser.add_argument("--quad-step", type=int, default=1)
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
    parser.add_argument("--q-values", default="1 2 3")
    parser.add_argument("--thresholds", default="1e-6 1e-5 1e-4 1e-3")
    parser.add_argument("--json-out", default="spectral_tail_moment_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = []
    print(
        f"Spectral tail moment certificate model={args.model} kind={args.kind}",
        flush=True,
    )
    print("  basis cons null gamma2      lambda_min  threshold tail_frac q1_bound q2_bound", flush=True)
    for basis in parse_ints(args.basis):
        row = compute_case(args, basis)
        if row is None:
            continue
        rows.append(serial(row))
        display_threshold = row["thresholds"][min(2, len(row["thresholds"]) - 1)]
        b1 = display_threshold["bounds"][0]["boundFrac"]
        b2 = display_threshold["bounds"][1]["boundFrac"] if len(display_threshold["bounds"]) > 1 else mp.mpf("0")
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['nullity']:4d} "
            f"{fmt(row['gamma2']):>10} {fmt(row['lambdaMin']):>10} "
            f"{fmt(display_threshold['threshold']):>9} {fmt(display_threshold['tailFrac']):>9} "
            f"{fmt(b1):>9} {fmt(b2):>9}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "qValues": [f(x) for x in parse_mpf_list(args.q_values)],
        "thresholds": [f(x) for x in parse_mpf_list(args.thresholds)],
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
