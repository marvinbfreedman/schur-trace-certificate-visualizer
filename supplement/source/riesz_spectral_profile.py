#!/usr/bin/env python3
r"""Spectral profile of the Riesz solve A w = h on ker(R).

For the quotient split V = ker(R) \oplus U, the top Douglas vector u gives

  h = B u,          w = A^+ h.

Here A is the form restricted to ker(R).  This script diagonalizes A and shows
which eigenmodes amplify the source into the endpoint boundary-layer witness.
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
    poly_value,
    quotient_certificate,
    trace_matrix,
)


def f(x):
    return float(x)


def fmt(x, digits=10):
    return mp.nstr(x, digits)


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


def split_spaces(K, R, rank_tol_text):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    return columns(rvecs, n_idx), columns(rvecs, u_idx), rvals, rank_tol


def compute(args):
    local = make_args(args)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, mp.mpf(local.omega), length)
    centers, R = trace_matrix(polys, local)
    cert = quotient_certificate(K, R, local)
    N, U, rvals, rank_tol = split_spaces(K, R, local.rank_tol)

    A = N.T * K * N
    B = N.T * K * U
    avals, avecs = mp.eigsy((A + A.T) / 2, eigvals_only=False)

    # Top eigenvector of Gamma^*Gamma = B^* A^+ B.
    inv_diag = [mp.mpf("0") if val <= mp.mpf(local.psd_tol) else 1 / val for val in avals]
    Aplus = avecs * mp.diag(inv_diag) * avecs.T
    gamma2 = B.T * Aplus * B
    gvals, gvecs = mp.eigsy((gamma2 + gamma2.T) / 2, eigvals_only=False)
    top = len(gvals) - 1
    u_coords = mp.matrix([gvecs[i, top] for i in range(gvecs.rows)])
    h_coords = B * u_coords

    # A-eigenbasis coordinates.
    h_eig = avecs.T * h_coords
    w_eig = mp.matrix([
        mp.mpf("0") if avals[i] <= mp.mpf(local.psd_tol) else h_eig[i] / avals[i]
        for i in range(len(avals))
    ])
    energy_terms = [
        mp.mpf("0") if avals[i] <= mp.mpf(local.psd_tol) else (h_eig[i] ** 2) / avals[i]
        for i in range(len(avals))
    ]
    witness_terms = [w_eig[i] ** 2 for i in range(len(avals))]
    source_terms = [h_eig[i] ** 2 for i in range(len(avals))]
    energy_total = mp.fsum(energy_terms)
    witness_total = mp.fsum(witness_terms)
    source_total = mp.fsum(source_terms)

    rows = []
    for i in range(len(avals)):
        mode_coeffs = N * columns(avecs, [i])
        grid = []
        for j in range(args.points):
            x = length * j / (args.points - 1)
            grid.append(f(eval_function(polys, [mode_coeffs[k, 0] for k in range(mode_coeffs.rows)], x)))
        peak_idx = max(range(len(grid)), key=lambda idx: abs(grid[idx]))
        rows.append(
            {
                "mode": i,
                "lambda": f(avals[i]),
                "sourceCoeff": f(h_eig[i]),
                "witnessCoeff": f(w_eig[i]),
                "sourceFrac": f(source_terms[i] / source_total) if source_total else 0.0,
                "witnessFrac": f(witness_terms[i] / witness_total) if witness_total else 0.0,
                "energyFrac": f(energy_terms[i] / energy_total) if energy_total else 0.0,
                "modePeakS": f(length * peak_idx / (args.points - 1)),
                "modePeakValue": grid[peak_idx],
            }
        )

    top_by_energy = sorted(rows, key=lambda row: row["energyFrac"], reverse=True)
    top_by_witness = sorted(rows, key=lambda row: row["witnessFrac"], reverse=True)
    return {
        "model": local.model,
        "kind": local.kind,
        "omega": f(mp.mpf(local.omega)),
        "L": f(length),
        "basis": local.basis,
        "constraints": local.constraints,
        "rank": cert["rank"],
        "nullity": cert["nullity"],
        "kMin": f(cert["kmin"]),
        "kerMin": f(cert["amin"]),
        "gamma2Top": f(gvals[top]),
        "energyTotal": f(energy_total),
        "witnessTotal": f(witness_total),
        "sourceTotal": f(source_total),
        "rows": rows,
        "topEnergyModes": top_by_energy[: args.top],
        "topWitnessModes": top_by_witness[: args.top],
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
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--json-out", default="riesz_spectral_profile.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    data = compute(args)
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")

    print(
        f"Riesz spectral profile model={data['model']} basis={data['basis']} "
        f"gamma2={fmt(data['gamma2Top'])}"
    )
    print("  top modes by energy contribution:")
    print("    mode  lambda        source      witness     energy_frac witness_frac peak_s")
    for row in data["topEnergyModes"]:
        print(
            f"    {row['mode']:4d} {fmt(row['lambda'], 8):>12} "
            f"{fmt(row['sourceCoeff'], 8):>11} {fmt(row['witnessCoeff'], 8):>11} "
            f"{row['energyFrac']:11.4f} {row['witnessFrac']:12.4f} "
            f"{row['modePeakS']:6.3f}"
        )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
