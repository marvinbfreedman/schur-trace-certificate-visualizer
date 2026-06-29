#!/usr/bin/env python3
r"""Finite candidate for the commuted Sturm elliptic estimate.

Endpoint energy A=K|ker(R) alone does not control high Sobolev graph norms.
The next natural commuted candidate is

    S_m = sum_{r=0}^m (D_s^r)^* K (D_s^r),

where K is the endpoint B kernel matrix on the full Galerkin basis and D_s is
the exact polynomial derivative operator in that basis.  This is the finite
analogue of commuting spectral derivatives through the endpoint Sturm Green
identity.

The script tests

    ||f||_{W^{m,2}}^2 <= C_{m,M} <S_m f,f>,
        f in H_M cap ker(R),

where H_M is defined by the low modes of A=K|ker(R).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import poly_derivative_coeffs, poly_inner, sobolev_matrix  # noqa: E402
from lagrange_split_control_certificate import parse_ints  # noqa: E402
from quotient_factorization_mp import columns, gram_matrix, trace_matrix  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def derivative_matrix(polys, length):
    """Matrix D with coefficients of d/ds in the orthonormal polynomial basis."""
    n = len(polys)
    D = mp.matrix(n)
    dpolys = [poly_derivative_coeffs(poly, 1) for poly in polys]
    for i, pi in enumerate(polys):
        for j, dpj in enumerate(dpolys):
            D[i, j] = poly_inner(pi, dpj, length)
    return D


def commuted_kernel_matrix(K, D, order: int):
    S = mp.matrix(K.rows)
    Dr = mp.eye(K.rows)
    for _r in range(order + 1):
        S += Dr.T * K * Dr
        Dr = D * Dr
    return (S + S.T) / 2


def positive_inverse_constant(numer, denom, tol):
    vals, vecs = mp.eigsy((denom + denom.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > tol]
    if not keep:
        return mp.inf, vals
    vkeep = columns(vecs, keep)
    scaled = mp.matrix(len(keep))
    denom_vals = [vals[i] for i in keep]
    core = vkeep.T * numer * vkeep
    for i in range(len(keep)):
        for j in range(len(keep)):
            scaled[i, j] = core[i, j] / mp.sqrt(denom_vals[i] * denom_vals[j])
    out_vals = mp.eigsy((scaled + scaled.T) / 2, eigvals_only=True)
    return max(mp.mpf("0"), out_vals[-1]), vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=11)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--cutoffs", default="2 3 4 5")
    parser.add_argument("--orders", default="2 4 6 8 10")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--denom-tol", default="1e-40")
    parser.add_argument("--json-out", default="lagrange_commuted_kernel_energy.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    qargs = make_qargs(args)
    K, polys = gram_matrix(qargs, mp.mpf(args.omega), mp.mpf(args.L))
    _centers, R = trace_matrix(polys, qargs)
    N, _U, rank, nullity = split_spaces(R, args.rank_tol)
    A = N.T * K * N
    avals_all, avecs_all = mp.eigsy((A + A.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(avals_all) if val > mp.mpf(args.psd_tol)]
    avals = [avals_all[i] for i in keep]
    a_modes = columns(avecs_all, keep)

    D = derivative_matrix(polys, mp.mpf(args.L))
    rows = []
    print(
        f"Commuted kernel energy basis={args.basis} constraints={args.constraints} "
        f"positive={len(avals)}",
        flush=True,
    )
    print("  cutoff order      C(W/S)      min(S)       max(S)", flush=True)
    for order in parse_ints(args.orders):
        W = N.T * sobolev_matrix(polys, mp.mpf(args.L), order, top_only=False) * N
        Sfull = commuted_kernel_matrix(K, D, order)
        S = N.T * Sfull * N
        for cutoff in parse_ints(args.cutoffs):
            cutoff = min(cutoff, len(avals))
            high = columns(a_modes, list(range(cutoff, len(avals))))
            Wh = high.T * W * high
            Sh = high.T * S * high
            const, svals = positive_inverse_constant(Wh, Sh, mp.mpf(args.denom_tol))
            rows.append(
                {
                    "cutoff": cutoff,
                    "order": order,
                    "constant": f(const),
                    "norm": f(mp.sqrt(max(mp.mpf("0"), const))) if const != mp.inf else None,
                    "sMin": f(svals[0]) if len(svals) else None,
                    "sMax": f(svals[-1]) if len(svals) else None,
                }
            )
            print(
                f"  {cutoff:6d} {order:5d} {fmt(const, 10):>12} "
                f"{fmt(svals[0] if len(svals) else 0, 8):>11} "
                f"{fmt(svals[-1] if len(svals) else 0, 8):>11}",
                flush=True,
            )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "constraints": args.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "lambdaMinA": f(avals[0]),
        "lambdaMaxA": f(avals[-1]),
        "rows": rows,
        "interpretation": (
            "S_m=sum_{r<=m}(D^r)^*K(D^r) is the finite commuted-kernel "
            "candidate for the Sturm elliptic estimate."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
