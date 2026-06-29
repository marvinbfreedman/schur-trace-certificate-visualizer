#!/usr/bin/env python3
r"""Closed-trace Hardy/Green representer certificate.

The continuum Hardy/Green inequality for the source-window residuals is

    ||E_u f||^2 + ||d_u E_u f||^2 <= C <A f,f>

on the closed-trace high block.  In Hilbert-space language this is a range
statement: every scalar component of E_u and d_u E_u must have an A-Green
representer g satisfying

    ell(f) = <g,f>_A,        ||ell||_{A^{-1}}^2 = ||g||_A^2.

This script constructs the finite Galerkin representers on H_M cap ker R,
checks the range identity in the A-eigenbasis, and records profile data for
the representers.  It is the finite model of the analytic Hardy/Green proof:
replace the matrix representers by solutions of the continuum Green problem
with the closed-trace boundary condition.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from aux_regularizer_certificate import direct_absorption_constant  # noqa: E402
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from lagrange_hardy_graph_certificate import residual_rows_for  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    poly_derivative_value,
    poly_value,
    trace_matrix,
)
from source_window_derivative_scan import residual_rows_du_for  # noqa: E402
from trace_lagrange_adjoint_control import load_exact_trace  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_floats(text: str) -> list[mp.mpf]:
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def row_component(row, comp: int):
    out = mp.matrix(1, row.cols)
    for j in range(row.cols):
        out[0, j] = row[comp, j]
    return out


def energy_representer(row_n_comp, A, high_a_modes, high_avals):
    coeffs = row_n_comp * high_a_modes
    rep_n = mp.matrix(high_a_modes.rows, 1)
    norm2 = mp.mpf("0")
    for j, lam in enumerate(high_avals):
        c = coeffs[0, j]
        norm2 += c * c / lam
        for i in range(high_a_modes.rows):
            rep_n[i] += high_a_modes[i, j] * c / lam

    predicted = rep_n.T * A * high_a_modes
    defect2 = mp.mpf("0")
    target2 = mp.mpf("0")
    for j in range(high_a_modes.cols):
        defect = predicted[0, j] - coeffs[0, j]
        defect2 += defect * defect
        target2 += coeffs[0, j] * coeffs[0, j]
    rel = mp.sqrt(defect2 / target2) if target2 else mp.mpf("0")
    return rep_n, norm2, rel


def full_coeffs(N, rep_n):
    out = N * rep_n
    return [out[i] for i in range(out.rows)]


def profile(polys, coeffs, length, points: int, max_deriv: int):
    rows = []
    max_abs = mp.mpf("0")
    max_derivs = [mp.mpf("0") for _ in range(max_deriv + 1)]
    l2 = mp.mpf("0")
    prev_x = None
    prev_val = None
    for i in range(points):
        x = length * i / (points - 1)
        derivs = []
        for d in range(max_deriv + 1):
            value = mp.fsum(
                coeffs[j] * poly_derivative_value(polys[j], x, d)
                for j in range(len(polys))
            )
            derivs.append(value)
            max_derivs[d] = max(max_derivs[d], abs(value))
        value = derivs[0]
        max_abs = max(max_abs, abs(value))
        if prev_x is not None:
            dx = x - prev_x
            l2 += dx * (prev_val * prev_val + value * value) / 2
        prev_x = x
        prev_val = value
        rows.append(
            {
                "s": f(x),
                "value": f(value),
                "d1": f(derivs[1]) if max_deriv >= 1 else 0.0,
                "d2": f(derivs[2]) if max_deriv >= 2 else 0.0,
            }
        )
    return {
        "maxAbs": f(max_abs),
        "l2Approx": f(l2),
        "maxDerivatives": [f(x) for x in max_derivs],
        "grid": rows,
    }


def component_records(args, polys, N, A, high_a_modes, high_avals, rows, label):
    records = []
    row_n = rows * N
    combined, _size = direct_absorption_constant(row_n, high_a_modes, high_avals)
    for comp, name in enumerate(["boundary", "adjointEval"]):
        row_comp = row_component(row_n, comp)
        rep_n, norm2, rel = energy_representer(row_comp, A, high_a_modes, high_avals)
        coeffs = full_coeffs(N, rep_n)
        prof = profile(
            polys,
            coeffs,
            mp.mpf(args.L),
            args.profile_points,
            args.profile_derivs,
        )
        records.append(
            {
                "label": label,
                "component": name,
                "norm2": f(norm2),
                "norm": f(mp.sqrt(max(mp.mpf("0"), norm2))),
                "rangeRelativeDefect": f(rel),
                "profile": prof,
            }
        )
    return f(combined), records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-values", default="0.08 0.30 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=20)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=80)
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
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--profile-points", type=int, default=81)
    parser.add_argument("--profile-derivs", type=int, default=2)
    parser.add_argument("--json-out", default="closed_trace_hardy_green_certificate.json")
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
    cutoff = min(args.cutoff, len(avals))
    high_a_modes = columns(a_modes, list(range(cutoff, len(avals))))
    high_avals = avals[cutoff:]

    _vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    rows_out = []
    records = []
    max_range_defect = mp.mpf("0")
    print(
        f"Closed-trace Hardy/Green certificate basis={args.basis} "
        f"constraints={args.constraints} cutoff={cutoff}",
        flush=True,
    )
    print("  u       E_const      dE_const     max_range_defect", flush=True)

    for u in parse_floats(args.source_values):
        e_rows, _pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, u)
        de_rows, _dpstar = residual_rows_du_for(args, polys, e_derivs, r_nodes, r_weights, u)
        e_const, e_records = component_records(
            args, polys, N, A, high_a_modes, high_avals, e_rows, "E"
        )
        de_const, de_records = component_records(
            args, polys, N, A, high_a_modes, high_avals, de_rows, "dE"
        )
        local_records = e_records + de_records
        for item in local_records:
            max_range_defect = max(max_range_defect, mp.mpf(item["rangeRelativeDefect"]))
        records.extend({"u": f(u), **item} for item in local_records)
        rows_out.append(
            {
                "u": f(u),
                "EConstant": e_const,
                "dEConstant": de_const,
                "combinedSumConstant": f(mp.mpf(e_const) + mp.mpf(de_const)),
                "maxRangeRelativeDefect": f(
                    max(mp.mpf(item["rangeRelativeDefect"]) for item in local_records)
                ),
            }
        )
        print(
            f"  {fmt(u, 6):>6} {fmt(mp.mpf(e_const), 10):>12} "
            f"{fmt(mp.mpf(de_const), 10):>12} {fmt(max_range_defect, 8):>16}",
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
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinHigh": f(high_avals[0]) if high_avals else None,
        "lambdaMaxHigh": f(high_avals[-1]) if high_avals else None,
        "maxRangeRelativeDefect": f(max_range_defect),
        "rows": rows_out,
        "representers": records,
        "theoremTemplate": (
            "If the continuum Green representers g_{u,k} for the components "
            "of E_u and partial_u E_u exist in the closed-trace high-block "
            "A-space with uniformly bounded A-norm, then the closed-trace "
            "Hardy/Green inequality follows by Cauchy-Schwarz."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
