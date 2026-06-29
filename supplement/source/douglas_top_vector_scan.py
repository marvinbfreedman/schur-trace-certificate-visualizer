#!/usr/bin/env python3
r"""Scan stability of the top Douglas eigenvector.

The top eigenvector of B^* A^+ B is only analytically useful if its shape is
stable as the Galerkin basis is refined.  This script computes the top vector
for several basis orders and reports simple shape diagnostics:

  - gamma2_top = lambda_max(B^* A^+ B)
  - peak location of the U-direction
  - peak location of the paired ker(R) witness A^+ B u
  - correlation with the previous normalized U-direction

It also writes a JSON file for the visual dashboard.
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

from douglas_top_vector import compute  # noqa: E402


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def constraint_count(args, basis):
    if args.constraint_rule == "half":
        return max(2, basis // 2)
    if args.constraint_rule == "basis":
        return basis
    return args.constraints


def local_args(args, basis):
    return SimpleNamespace(
        model=args.model,
        kind=args.kind,
        omega=args.omega,
        L=args.L,
        basis=basis,
        quad=args.quad_base + args.quad_step * basis,
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=constraint_count(args, basis),
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
        points=args.points,
        json_out="",
    )


def max_abs_point(grid, key):
    idx = max(range(len(grid)), key=lambda i: abs(grid[i][key]))
    return grid[idx]["s"], grid[idx][key]


def sign_changes(values):
    count = 0
    prev = 0
    for value in values:
        sign = 1 if value > 0 else -1 if value < 0 else 0
        if sign and prev and sign != prev:
            count += 1
        if sign:
            prev = sign
    return count


def normalized(values):
    scale = max(max(abs(v) for v in values), 1e-300)
    return [v / scale for v in values]


def corr(a, b):
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot / (na * nb)


def summarize(data, prev_u):
    u_values = [row["u"] for row in data["grid"]]
    n_values = [row["nWitness"] for row in data["grid"]]
    u_norm = normalized(u_values)
    n_norm = normalized(n_values)
    if prev_u is not None and corr(prev_u, u_norm) < 0:
        u_norm = [-x for x in u_norm]
        n_norm = [-x for x in n_norm]
    u_peak_s, u_peak_v = max_abs_point(data["grid"], "u")
    n_peak_s, n_peak_v = max_abs_point(data["grid"], "nWitness")
    trace_abs = max(abs(x) for x in data["traceValues"])
    return {
        "basis": data["basis"],
        "constraints": data["constraints"],
        "quad": data["quad"],
        "gamma2Top": data["gamma2Top"],
        "gammaTop": data["gammaTop"],
        "kMin": data["kMin"],
        "kerMin": data["kerMin"],
        "rangeResidual": data["rangeResidual"],
        "uPeakS": u_peak_s,
        "uPeakValue": u_peak_v,
        "nPeakS": n_peak_s,
        "nPeakValue": n_peak_v,
        "uSignChanges": sign_changes(u_values),
        "nSignChanges": sign_changes(n_values),
        "maxTraceAbs": trace_abs,
        "corrPrev": None if prev_u is None else corr(prev_u, u_norm),
        "uNorm": u_norm,
        "nNorm": n_norm,
        "gridS": [row["s"] for row in data["grid"]],
    }


def fmt(x):
    if x is None:
        return "n/a"
    if abs(x) >= 1e-3 and abs(x) < 1e4:
        return f"{x:.6g}"
    return f"{x:.3e}"


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
    parser.add_argument("--points", type=int, default=81)
    parser.add_argument("--json-out", default="douglas_top_vector_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    basis_values = parse_ints(args.basis)
    rows = []
    prev_u = None
    print(
        f"top Douglas vector scan model={args.model} kind={args.kind} "
        f"omega={args.omega} L={args.L}",
        flush=True,
    )
    print(
        "  basis cons gamma2_top  corr_prev  u_peak_s  n_peak_s  "
        "u_zeros n_zeros max|Ru|",
        flush=True,
    )
    for basis in basis_values:
        data = compute(local_args(args, basis))
        row = summarize(data, prev_u)
        rows.append(row)
        prev_u = row["uNorm"]
        print(
            f"  {row['basis']:5d} {row['constraints']:4d} "
            f"{fmt(row['gamma2Top']):>10} {fmt(row['corrPrev']):>10} "
            f"{row['uPeakS']:8.4f} {row['nPeakS']:8.4f} "
            f"{row['uSignChanges']:7d} {row['nSignChanges']:7d} "
            f"{fmt(row['maxTraceAbs']):>10}",
            flush=True,
        )

    out = {
        "model": args.model,
        "kind": args.kind,
        "omega": float(args.omega),
        "L": float(args.L),
        "basisValues": basis_values,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
