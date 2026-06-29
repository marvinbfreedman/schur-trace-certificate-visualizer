#!/usr/bin/env python3
"""Scan finite Douglas constants under Galerkin refinement.

The continuum target is the uniform estimate

  |b(n,u)|^2 <= C_D a(n,n) ||u||_U^2,

which is measured in finite dimension by

  ||Gamma||^2 = lambda_max(B^* A^+ B).

This script reuses quotient_factorization_mp.py and runs a short refinement
table.  It is a diagnostic, not a proof: growth in the reported constants
would point to failure or a missing norm; stability points to the analytic
inequality that should be proved.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from quotient_factorization_mp import (  # noqa: E402
    gram_matrix,
    quotient_certificate,
    trace_matrix,
)


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def fmt(x, digits=10):
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="10 12 14 16")
    parser.add_argument("--constraint-rule", choices=["half", "fixed", "basis"], default="half")
    parser.add_argument("--constraints", type=int, default=8)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--quad-base", type=int, default=6)
    parser.add_argument("--quad-step", type=int, default=1)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=14)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=22)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=55)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    basis_values = parse_ints(args.basis)

    print(
        f"Douglas refinement scan model={args.model} kind={args.kind} "
        f"omega={args.omega} L={args.L} dps={args.dps}",
        flush=True,
    )
    print(
        "  basis  cons  quad  rank null  K_min         ker_min       "
        "range_resid    gamma2_max    H_max         repaired_min",
        flush=True,
    )

    rows = []
    for basis in basis_values:
        if args.constraint_rule == "half":
            constraints = max(2, basis // 2)
        elif args.constraint_rule == "basis":
            constraints = basis
        else:
            constraints = args.constraints
        local = make_args(args, basis, constraints)
        K, polys = gram_matrix(local, mp.mpf(local.omega), mp.mpf(local.L))
        _, R = trace_matrix(polys, local)
        cert = quotient_certificate(K, R, local)
        rows.append(cert)
        print(
            f"  {basis:5d} {constraints:5d} {local.quad:5d} "
            f"{cert['rank']:5d} {cert['nullity']:4d} "
            f"{fmt(cert['kmin'], 8):>12} {fmt(cert['amin'], 8):>12} "
            f"{fmt(cert['normalized_range_resid'], 8):>12} "
            f"{fmt(cert['gamma2_max'], 8):>12} "
            f"{fmt(cert['h_max'], 8):>12} {fmt(cert['p2min'], 8):>12}",
            flush=True,
        )

    gamma_max = max(row["gamma2_max"] for row in rows)
    gamma_min = min(row["gamma2_max"] for row in rows)
    print(f"  gamma2 range: min={fmt(gamma_min, 12)} max={fmt(gamma_max, 12)}", flush=True)


if __name__ == "__main__":
    main()
