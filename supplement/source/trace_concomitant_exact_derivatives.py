#!/usr/bin/env python3
r"""Use confluent-kernel identities to differentiate the trace row analytically.

The finite-difference trace tests suggested that the Green concomitant of

    Lambda_s(f)=sum_{k=0}^8 e_k(s) f^(k)(s)/k!

almost lies in the span of ``D_s^q Lambda_s``.  This script removes the
finite-difference step.  It obtains derivatives of the endpoint defect row
``e(s)`` from the confluent eigenvalue equation

    J(s)e(s)=lambda(s)e(s),

where ``J(s)`` is the Taylor-normalized confluent matrix of the endpoint
kernel.  The center derivatives of ``J`` are exact consequences of the larger
confluent matrix:

    [z^i w^j delta^p] K(s+delta+z,s+delta+w)
      = sum_{a+b=p} binom(i+a,i) binom(j+b,j) J_{i+a,j+b}(s).

The script then rebuilds the trace Green concomitant row and repeats the
membership test.  It also prints the algebraic obstruction: since the leading
coefficient ``e_8(s)/8!`` is nonzero, no nonzero row of order at most 7 can be
exactly equal to ``M Lambda`` for an ordinary differential operator ``M``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_kb_confluent_mp import integrate  # noqa: E402
from lambda_differential_closure import trace_derivative_row  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_concomitant_membership import solve_membership, source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def vec_col(mat, col: int) -> list[mp.mpf]:
    return [mat[i, col] for i in range(mat.rows)]


def dot(a, b):
    return mp.fsum(x * y for x, y in zip(a, b))


def mat_vec(mat, vec):
    return [mp.fsum(mat[i, j] * vec[j] for j in range(mat.cols)) for i in range(mat.rows)]


def vec_add(a, b):
    return [x + y for x, y in zip(a, b)]


def vec_scale(a, c):
    return [c * x for x in a]


def center_taylor_mats(big, base_order: int, max_q: int):
    mats = []
    for p in range(max_q + 1):
        mat = mp.matrix(base_order)
        for i in range(base_order):
            for j in range(base_order):
                mat[i, j] = mp.fsum(
                    mp.binomial(i + a, i)
                    * mp.binomial(j + p - a, j)
                    * big[i + a, j + p - a]
                    for a in range(p + 1)
                )
        mats.append((mat + mat.T) / 2)
    return mats


def eigen_taylor(A_coeffs, max_q: int):
    vals, vecs = mp.eigsy(A_coeffs[0], eigvals_only=False)
    idx = 0
    v0 = vec_col(vecs, idx)
    if v0[0] > 0:
        v0 = [-x for x in v0]
    eigvecs = [vec_col(vecs, j) for j in range(vecs.cols)]

    lam = [vals[idx]]
    vec = [v0]
    for n in range(1, max_q + 1):
        rhs = [mp.mpf("0") for _ in v0]
        for p in range(1, n + 1):
            rhs = vec_add(rhs, vec_scale(mat_vec(A_coeffs[p], vec[n - p]), -1))
        for p in range(1, n):
            rhs = vec_add(rhs, vec_scale(vec[n - p], lam[p]))

        lam_n = -dot(v0, rhs)
        rhs = vec_add(rhs, vec_scale(v0, lam_n))

        part = [mp.mpf("0") for _ in v0]
        for j, qj in enumerate(eigvecs):
            if j == idx:
                continue
            coeff = dot(qj, rhs) / (vals[j] - vals[idx])
            part = vec_add(part, vec_scale(qj, coeff))

        beta = -mp.mpf("0.5") * mp.fsum(dot(vec[p], vec[n - p]) for p in range(1, n))
        vec_n = vec_add(part, vec_scale(v0, beta))
        vec.append(vec_n)
        lam.append(lam_n)

    e_derivs = [
        [mp.factorial(q) * value for value in vec[q]]
        for q in range(max_q + 1)
    ]
    lam_derivs = [mp.factorial(q) * value for q, value in enumerate(lam)]
    return vals, e_derivs, lam_derivs


def pad_row(row, cols):
    return row + [mp.mpf("0") for _ in range(cols - len(row))]


def obstruction_text(e_derivs, h0, target):
    leading = e_derivs[0][-1] / mp.factorial(len(e_derivs[0]) - 1)
    predicted = -leading * h0
    return {
        "a8": leading,
        "h0": h0,
        "targetF7": target[-1],
        "predictedF7": predicted,
        "defect": target[-1] - predicted,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--max-trace-q", type=int, default=10)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--json-out", default="trace_concomitant_exact_derivatives.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    big_order = args.jet_order + args.max_trace_q

    big, _segments = integrate("kb", c, s0, big_order, mp.mpf(args.matrix_rmax), args.matrix_order)
    A_coeffs = center_taylor_mats(big, args.jet_order, args.max_trace_q)
    vals, e_derivs, lam_derivs = eigen_taylor(A_coeffs, args.max_trace_q)
    trace_rows = [
        trace_derivative_row(e_derivs, q, args.jet_order)
        for q in range(args.max_trace_q + 1)
    ]

    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    rows = []
    obstructions = []

    print(
        f"Trace concomitant exact derivatives s0={s0} "
        f"max_trace_q={args.max_trace_q} big_order={big_order}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} "
        f"gap={fmt(vals[1]-vals[0], 12)} min|e8|={fmt(abs(e_derivs[0][-1]), 12)}",
        flush=True,
    )
    print("  t      qmax  generators  rel_residual  alpha_norm  target_norm", flush=True)

    for t in t_values:
        h_derivs = source_derivatives(
            c,
            s0,
            t,
            args.jet_order - 2,
            r_nodes,
            r_weights,
        )
        target = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
        obstructions.append({"t": f(t), **{k: f(v) for k, v in obstruction_text(e_derivs, h_derivs[0], target).items()}})
        for qmax in range(args.max_trace_q + 1):
            generators = trace_rows[: qmax + 1]
            cols = max(len(target), max(len(row) for row in generators))
            alpha, res_norm, target_norm, rel = solve_membership(
                pad_row(target, cols),
                [pad_row(row, cols) for row in generators],
            )
            alpha_norm = mp.sqrt(mp.fsum(a * a for a in alpha)) if alpha else mp.mpf("0")
            rows.append(
                {
                    "t": f(t),
                    "maxTraceQ": qmax,
                    "generators": len(generators),
                    "relativeResidual": f(rel),
                    "residualNorm": f(res_norm),
                    "targetNorm": f(target_norm),
                    "alphaNorm": f(alpha_norm),
                    "alphaMaxAbs": f(max(abs(a) for a in alpha)) if alpha else 0.0,
                }
            )
            if qmax in (0, 6, 8, 9, 10):
                print(
                    f"  {fmt(t, 6):>6} {qmax:5d} {len(generators):10d} "
                    f"{fmt(rel, 8):>12} {fmt(alpha_norm, 8):>11} {fmt(target_norm, 8):>11}",
                    flush=True,
                )

    data = {
        "s0": f(s0),
        "jetOrder": args.jet_order,
        "maxTraceQ": args.max_trace_q,
        "bigOrder": big_order,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "obstructions": obstructions,
        "algebraicConclusion": (
            "Exact left-ideal membership is impossible for a nonzero order<=7 row "
            "because e_8(s)/8! is nonzero."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("  obstruction coefficient:")
    for item in obstructions:
        print(
            f"    t={item['t']} a8={fmt(item['a8'], 8)} h0={fmt(item['h0'], 8)} "
            f"target_f7={fmt(item['targetF7'], 8)}",
            flush=True,
        )
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
