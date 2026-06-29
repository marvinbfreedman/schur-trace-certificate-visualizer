#!/usr/bin/env python3
r"""Finite spectral diagnostic for the commuted source estimate.

The abstract high-tail reduction needs a bound of the form

    || A^{m/2} B u ||^2 <= M_m^2 C_full ||u||^2.

In the finite quotient model, A=K|ker(R) and B is the cross block from trace
directions into ker(R).  If A e_j=lambda_j e_j and b_j(u)=<B u,e_j>, this
script computes

    || A^{m/2} B ||^2 = || diag(lambda_j^{m/2}) B_eig ||^2

and the Hilbert-Schmidt scalar moment

    sum_j lambda_j^m ||b_j||^2.

The numbers do not prove the continuum commuted estimate; they check whether
the proposed integration-by-parts gain is compatible with the same Galerkin
sections used in the source-envelope certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from high_frequency_tail_refinement import constraints_for  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, max_eig_or_zero, trace_matrix  # noqa: E402


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def make_args(base, basis: int, constraints: int):
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


def split_spaces(R, rank_tol_text: str):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx), len(u_idx), len(n_idx)


def row_norm2(mat, row_index: int):
    return mp.fsum(mat[row_index, j] ** 2 for j in range(mat.cols))


def weighted_source_norm2(beig, avals, moment: float):
    if beig.rows == 0 or beig.cols == 0:
        return mp.mpf("0")
    mat = mp.matrix(beig.cols)
    for j, lam in enumerate(avals):
        weight = lam ** mp.mpf(moment)
        for a in range(beig.cols):
            for b in range(beig.cols):
                mat[a, b] += weight * beig[j, a] * beig[j, b]
    return max_eig_or_zero(mat)


def weighted_source_hs2(beig, avals, moment: float):
    return mp.fsum((lam ** mp.mpf(moment)) * row_norm2(beig, j) for j, lam in enumerate(avals))


def tail_start_for_lambda(avals, threshold):
    for j, lam in enumerate(avals):
        if lam >= threshold:
            return j
    return len(avals)


def window_constant(beig, avals, start: int, stop: int, inverse_power=1):
    if beig.cols == 0 or start >= stop:
        return mp.mpf("0")
    mat = mp.matrix(beig.cols)
    for j in range(start, stop):
        scale = avals[j] ** (-mp.mpf(inverse_power))
        for a in range(beig.cols):
            for b in range(beig.cols):
                mat[a, b] += scale * beig[j, a] * beig[j, b]
    return max_eig_or_zero(mat)


def compute_case(base, basis: int, moments: list[float], thresholds: list[float]):
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
    gamma2 = window_constant(beig, avals, 0, len(avals), inverse_power=1)

    moment_rows = []
    for moment in moments:
        op = weighted_source_norm2(beig, avals, moment)
        hs = weighted_source_hs2(beig, avals, moment)
        tail_bounds = []
        for threshold in thresholds:
            lam0 = mp.mpf(threshold)
            start = tail_start_for_lambda(avals, lam0)
            actual = window_constant(beig, avals, start, len(avals), inverse_power=1)
            # From ||A^{m/2}B||^2: C_tail(Lambda)<=Lambda^{-(m+1)} ||A^{m/2}B||^2.
            op_bound = (lam0 ** (-(mp.mpf(moment) + 1))) * op
            hs_bound = (lam0 ** (-(mp.mpf(moment) + 1))) * hs
            tail_bounds.append(
                {
                    "lambdaThreshold": threshold,
                    "start": start,
                    "actualTailFrac": actual / gamma2 if gamma2 else mp.mpf("0"),
                    "operatorBoundFrac": op_bound / gamma2 if gamma2 else mp.mpf("0"),
                    "hsBoundFrac": hs_bound / gamma2 if gamma2 else mp.mpf("0"),
                }
            )
        moment_rows.append(
            {
                "moment": moment,
                "operatorNorm2": op,
                "operatorFrac": op / gamma2 if gamma2 else mp.mpf("0"),
                "hilbertSchmidt2": hs,
                "hilbertSchmidtFrac": hs / gamma2 if gamma2 else mp.mpf("0"),
                "tailBounds": tail_bounds,
            }
        )

    return {
        "basis": basis,
        "constraints": constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "gamma2": gamma2,
        "lambdaMin": avals[0],
        "lambdaMax": avals[-1],
        "moments": moment_rows,
    }


def serial(row: dict):
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
                "moment": item["moment"],
                "operatorNorm2": f(item["operatorNorm2"]),
                "operatorFrac": f(item["operatorFrac"]),
                "hilbertSchmidt2": f(item["hilbertSchmidt2"]),
                "hilbertSchmidtFrac": f(item["hilbertSchmidtFrac"]),
                "tailBounds": [
                    {
                        "lambdaThreshold": tb["lambdaThreshold"],
                        "start": tb["start"],
                        "actualTailFrac": f(tb["actualTailFrac"]),
                        "operatorBoundFrac": f(tb["operatorBoundFrac"]),
                        "hsBoundFrac": f(tb["hsBoundFrac"]),
                    }
                    for tb in item["tailBounds"]
                ],
            }
            for item in row["moments"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="12 14 16 18 20")
    parser.add_argument("--moments", default="0 0.5 1 2")
    parser.add_argument("--thresholds", default="1e-4 1e-3 1e-2 1e-1")
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
    parser.add_argument("--json-out", default="commuted_source_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    moments = parse_floats(args.moments)
    thresholds = parse_floats(args.thresholds)
    rows = []

    print(
        f"Commuted source certificate model={args.model} kind={args.kind} moments={moments}",
        flush=True,
    )
    print("  basis cons null Gamma^2     lambda range        m=1 op/G   m=1 hs/G", flush=True)
    for basis in parse_ints(args.basis):
        row = compute_case(args, basis, moments, thresholds)
        if row is None:
            continue
        rows.append(serial(row))
        m1 = min(row["moments"], key=lambda item: abs(item["moment"] - 1.0))
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['nullity']:4d} "
            f"{fmt(row['gamma2']):>10} "
            f"[{fmt(row['lambdaMin'], 4)}, {fmt(row['lambdaMax'], 4)}] "
            f"{fmt(m1['operatorFrac']):>10} {fmt(m1['hilbertSchmidtFrac']):>10}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "moments": moments,
        "thresholds": thresholds,
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
