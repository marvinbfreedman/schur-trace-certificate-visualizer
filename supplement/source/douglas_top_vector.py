#!/usr/bin/env python3
r"""Compute the top Douglas eigenvector.

For the quotient split V = ker(R) \oplus U, the finite Douglas operator is

  Gamma^* Gamma = B^* A^+ B,

where A is K restricted to ker(R) and B is the cross block.  The top
eigenvector of this operator is the U-direction that makes the Douglas bound

  |b(n,u)|^2 <= C_D a(n,n) ||u||^2

sharp in the finite model.

This script computes that vector, evaluates it as a function of s in the
shifted Legendre basis, and also evaluates the paired ker(R) witness

  n_u = A^+ B u.
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
    minmax_eigs,
    poly_value,
    positive_part_inverse,
    quotient_certificate,
    trace_matrix,
)


def fmt(x, digits=14):
    return mp.nstr(x, digits)


def f(x):
    return float(x)


def eval_function(polys, coeffs, x):
    return mp.fsum(coeffs[j] * poly_value(polys[j], x) for j in range(len(polys)))


def make_args(args):
    return SimpleNamespace(
        model=args.model,
        kind=args.kind,
        omega=args.omega,
        L=args.L,
        basis=args.basis,
        quad=args.quad,
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=args.constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
        margin=args.margin,
        dps=args.dps,
    )


def compute(args):
    local = make_args(args)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, mp.mpf(local.omega), length)
    centers, R = trace_matrix(polys, local)
    cert = quotient_certificate(K, R, local)

    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(local.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)

    A = N.T * K * N
    B = N.T * K * U
    _, _, _, _, Aplus = positive_part_inverse(A, mp.mpf(local.psd_tol))
    gamma2 = B.T * Aplus * B if B.cols else mp.matrix(0)
    gvals, gvecs = mp.eigsy((gamma2 + gamma2.T) / 2, eigvals_only=False)
    top = len(gvals) - 1
    u_coords = mp.matrix([gvecs[i, top] for i in range(gvecs.rows)])
    full_u = U * u_coords
    n_coords = Aplus * B * u_coords
    full_n = N * n_coords

    # Stable display sign: make the largest coefficient of u positive.
    max_idx = max(range(full_u.rows), key=lambda i: abs(full_u[i]))
    if full_u[max_idx] < 0:
        full_u = -full_u
        full_n = -full_n
        u_coords = -u_coords
        n_coords = -n_coords

    source_coords = B * u_coords
    full_source = N * source_coords
    trace_values = R * full_u
    grid = []
    points = args.points
    for i in range(points):
        x = length * i / (points - 1)
        grid.append(
            {
                "s": f(x),
                "u": f(eval_function(polys, full_u, x)),
                "nWitness": f(eval_function(polys, full_n, x)),
                "source": f(eval_function(polys, full_source, x)),
            }
        )

    return {
        "model": local.model,
        "kind": local.kind,
        "omega": f(mp.mpf(local.omega)),
        "L": f(length),
        "basis": local.basis,
        "constraints": local.constraints,
        "quad": local.quad,
        "endpointKernelOrder": local.endpoint_kernel_order,
        "endpointOrder": local.endpoint_order,
        "rank": cert["rank"],
        "nullity": cert["nullity"],
        "gamma2Top": f(gvals[top]),
        "gammaTop": f(mp.sqrt(max(mp.mpf("0"), gvals[top]))),
        "gamma2Eigenvalues": [f(x) for x in gvals],
        "kMin": f(cert["kmin"]),
        "kerMin": f(cert["amin"]),
        "rangeResidual": f(cert["normalized_range_resid"]),
        "uCoefficients": [f(full_u[i]) for i in range(full_u.rows)],
        "nWitnessCoefficients": [f(full_n[i]) for i in range(full_n.rows)],
        "sourceCoefficients": [f(full_source[i]) for i in range(full_source.rows)],
        "traceCenters": [f(x) for x in centers],
        "traceValues": [f(trace_values[i]) for i in range(trace_values.rows)],
        "grid": grid,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=16)
    parser.add_argument("--quad", type=int, default=22)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--constraints", type=int, default=8)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--points", type=int, default=121)
    parser.add_argument("--json-out", default="douglas_top_vector.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    data = compute(args)
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(
        f"top Douglas vector model={data['model']} basis={data['basis']} "
        f"constraints={data['constraints']}"
    )
    print(f"  gamma2_top={fmt(data['gamma2Top'])}")
    print(f"  gamma_top={fmt(data['gammaTop'])}")
    print(f"  K_min={fmt(data['kMin'])} K|kerR_min={fmt(data['kerMin'])}")
    print(f"  range_resid={fmt(data['rangeResidual'])}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
