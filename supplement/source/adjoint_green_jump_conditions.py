#!/usr/bin/env python3
r"""Derive the interior jump conditions for the actual adjoint Green kernel.

Let

    P f(a) = Lambda_a(f) = sum_{k=0}^8 a_k(a) f^(k)(a),
    a_k(a)=e_k(a)/k!,

and

    P^* K = sum_{k=0}^8 (-1)^k D_a^k(a_k K).

The source rows used in the proof are local jet functionals at s0:

    ell(f)=sum_{q=0}^7 d_q f^(q)(s0).

Therefore the actual adjoint Green coefficient K is not globally smooth.  If K
is C^8 away from s0 and has jumps

    Delta_r = K^(r)(s0+) - K^(r)(s0-),       0 <= r <= 7,

then the singular part of P^*K at s0 is

    sum_{q=0}^7 S_q delta_s0^(q),

where

    S_q = sum_{k=q+1}^8 (-1)^k
          [D^(k-1-q)(a_k K)]_{s0}.

To represent ell, one needs S_q=(-1)^q d_q.  This script constructs the exact
8x8 jump matrix Delta -> S and solves it for the two source components

    ell_1(f)=B_P[h_u,f](s0),
    ell_2(f)=(P^*h_u)(s0) f(s0).

This closes the distributional/source part of the adjoint Green derivation.
The still-open analytic part is the endpoint concomitant and compactness of the
piecewise homogeneous solutions in the trace-dual norm.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 10) -> str:
    return mp.nstr(x, digits)


def row_norm(values) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in values))


def coefficient_derivative(e_derivs, deriv: int, k: int) -> mp.mpf:
    return e_derivs[deriv][k] / mp.factorial(k)


def jump_matrix(e_derivs, order: int = 8) -> mp.matrix:
    mat = mp.matrix(order, order)
    for q in range(order):
        for r in range(order):
            total = mp.mpf("0")
            for k in range(q + 1, order + 1):
                n = k - 1 - q
                if r <= n:
                    total += (
                        ((-1) ** k)
                        * mp.binomial(n, r)
                        * coefficient_derivative(e_derivs, n - r, k)
                    )
            mat[q, r] = total
    return mat


def target_from_source_row(source_row):
    return [((-1) ** q) * source_row[q] for q in range(len(source_row))]


def solve_jump(mat: mp.matrix, source_row):
    target = target_from_source_row(source_row)
    rhs = mp.matrix(len(target), 1)
    for i, value in enumerate(target):
        rhs[i] = value
    delta = mp.lu_solve(mat, rhs)
    residual = mat * delta - rhs
    return (
        [delta[i] for i in range(delta.rows)],
        [residual[i] for i in range(residual.rows)],
        target,
    )


def mat_frobenius(mat: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=8)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--json-out", default="adjoint_green_jump_conditions.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, args.max_trace_q)
    jmat = jump_matrix(e_derivs, args.jet_order - 1)
    jinv = jmat ** -1
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    rows = []
    max_residual = mp.mpf("0")
    max_jump_norm = mp.mpf("0")
    for t in t_values:
        h_derivs = source_derivatives(
            c,
            s0,
            t,
            args.jet_order,
            r_nodes,
            r_weights,
        )
        boundary_row = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
        pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
        eval_row = [pstar] + [mp.mpf("0") for _ in range(args.jet_order - 2)]
        for label, source_row in (("boundary", boundary_row), ("adjointEval", eval_row)):
            delta, residual, target = solve_jump(jmat, source_row)
            res_norm = row_norm(residual)
            src_norm = row_norm(source_row)
            jump_norm = row_norm(delta)
            max_residual = max(max_residual, res_norm)
            max_jump_norm = max(max_jump_norm, jump_norm)
            rows.append(
                {
                    "sourceNode": f(t),
                    "component": label,
                    "sourceRowNorm": f(src_norm),
                    "targetNorm": f(row_norm(target)),
                    "jumpNorm": f(jump_norm),
                    "maxAbsJump": f(max(abs(x) for x in delta)),
                    "residualNorm": f(res_norm),
                    "relativeResidual": f(res_norm / max(mp.mpf("1"), row_norm(target))),
                    "jumpVector": [f(x) for x in delta],
                    "sourceRow": [f(x) for x in source_row],
                }
            )

    data = {
        "theoremName": "adjoint Green interior jump conditions",
        "s0": f(s0),
        "jetOrder": args.jet_order,
        "operatorOrder": args.jet_order - 1,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "jumpMatrixFrobenius": f(mat_frobenius(jmat)),
        "inverseJumpMatrixFrobenius": f(mat_frobenius(jinv)),
        "jumpMatrixDiagonalProduct": f(mp.det(jmat)),
        "maxJumpSolveResidual": f(max_residual),
        "maxJumpNorm": f(max_jump_norm),
        "jumpMatrix": [[f(jmat[i, j]) for j in range(jmat.cols)] for i in range(jmat.rows)],
        "rows": rows,
        "closedStatements": {
            "distributionalSourceJumpLaw": True,
            "jumpMatrixInvertible": abs(mp.det(jmat)) > mp.mpf("1e-80"),
            "interiorSourceSolved": max_residual < mp.mpf("1e-50"),
            "endpointConcomitantSolved": False,
        },
        "formula": {
            "sourceDistribution": "ell(f)=sum_q d_q f^(q)(s0) is represented by sum_q (-1)^q d_q delta_s0^(q)",
            "jumpCondition": "sum_{k=q+1}^8 (-1)^k [D^(k-1-q)(a_k K)]_{s0}=(-1)^q d_q",
            "piecewiseEquation": "P^*K=0 on (a_-,s0) union (s0,a_+) after the source jump is imposed",
            "remainingBoundaryProblem": "choose the homogeneous constants so the endpoint concomitant is killed or lies in the source-inactive tail, with uniform trace-dual bounds",
        },
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Adjoint Green interior jump conditions", flush=True)
    print(f"  det jump matrix: {fmt(mp.det(jmat), 12)}", flush=True)
    print(f"  ||J^-1||_F: {fmt(mat_frobenius(jinv), 12)}", flush=True)
    print(f"  max jump solve residual: {fmt(max_residual, 12)}", flush=True)
    print(f"  max jump norm: {fmt(max_jump_norm, 12)}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
