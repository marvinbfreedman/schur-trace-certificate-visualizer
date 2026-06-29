#!/usr/bin/env python3
r"""Show what the closed trace condition alone does and does not imply.

The continuum condition

    Lambda_s(f)=0 for s in an interval

is a variable-coefficient eighth-order ODE when e_8(s) != 0.  It recursively
determines higher jets from the eight initial jets f, ..., f^(7).  Therefore it
does not by itself imply f=0, nor does it by itself kill arbitrary local source
rows depending on those initial jets.

This script constructs the finite jet tower at s0 from exact trace derivatives
and projects the Lagrange source rows onto the local closed-trace jet solution
space.  Nonzero projection means the uniqueness lemma must use the Volterra/
Sturm energy equation or prove source-active range inclusion, not just the ODE
Lambda_s(f)=0.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from lagrange_hardy_graph_certificate import residual_rows_for  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature, shifted_legendre_polys  # noqa: E402
from trace_lagrange_adjoint_control import load_exact_trace  # noqa: E402
from lambda_differential_closure import trace_derivative_row  # noqa: E402


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_floats(text: str) -> list[mp.mpf]:
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def solve_closed_trace_jets(e_derivs, jet_order: int, max_q: int):
    """Return matrix mapping initial jets f0..f7 to jets f0..f_(7+max_q+1)."""
    initial = jet_order - 1
    cols = initial
    total_jets = initial + max_q + 1
    out = mp.matrix(total_jets, cols)
    for i in range(initial):
        out[i, i] = 1
    for q in range(max_q + 1):
        row = trace_derivative_row(e_derivs, q, jet_order)
        pivot = initial + q
        rhs_coeffs = [mp.mpf("0") for _ in range(cols)]
        for j in range(pivot):
            for c in range(cols):
                rhs_coeffs[c] -= row[j] * out[j, c]
        for c in range(cols):
            out[pivot, c] = rhs_coeffs[c] / row[pivot]
    return out


def jet_eval_rows(args, e_derivs, t):
    """Compute the two source rows on monomial Taylor jets at s0."""
    # Reuse residual_rows_for by feeding Taylor monomials p_j(s)=(s-s0)^j.
    # Then p_j^(k)(s0)=j! if j=k else 0, so divide by factorial(j) to get
    # rows acting on ordinary jet coordinates f^(j)(s0).
    s0 = mp.mpf(args.s0)
    degree = args.jet_order + args.max_trace_q
    polys = []
    for j in range(degree):
        coeffs = [mp.mpf("0") for _ in range(j + 1)]
        for k in range(j + 1):
            coeffs[k] = mp.binomial(j, k) * (-s0) ** (j - k)
        polys.append(coeffs)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    rows, pstar = residual_rows_for(args, polys, e_derivs, r_nodes, r_weights, t)
    out = mp.matrix(rows.rows, rows.cols)
    for i in range(rows.rows):
        for j in range(rows.cols):
            out[i, j] = rows[i, j] / mp.factorial(j)
    return out, pstar


def row_norm(mat):
    return mp.sqrt(mp.fsum(mat[i, j] ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--t-values", default="0.08 0.30 0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=12)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="closed_trace_local_ode_obstruction.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, args.max_trace_q)
    jet_map = solve_closed_trace_jets(e_derivs, args.jet_order, args.max_trace_q)
    t_values = parse_floats(args.t_values)

    rows = []
    print(
        f"Closed trace local ODE obstruction s0={args.s0} max_q={args.max_trace_q}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} lambda1={fmt(vals[1], 12)} e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  t       projected_source_norm  pstar       eval_projection_norm", flush=True)
    for t in t_values:
        source_rows, pstar = jet_eval_rows(args, e_derivs, t)
        projected = source_rows * jet_map
        eval_projection_norm = abs(projected[1, 0]) if projected.rows > 1 else mp.mpf("0")
        norm = row_norm(projected)
        rows.append(
            {
                "t": f(t),
                "pStar": f(pstar),
                "projectedSourceNorm": f(norm),
                "evalProjectionNorm": f(eval_projection_norm),
                "projectedRows": [[f(projected[i, j]) for j in range(projected.cols)] for i in range(projected.rows)],
            }
        )
        print(
            f"  {fmt(t, 6):>6} {fmt(norm, 10):>22} {fmt(pstar, 8):>11} {fmt(eval_projection_norm, 10):>22}",
            flush=True,
        )

    data = {
        "s0": f(mp.mpf(args.s0)),
        "jetOrder": args.jet_order,
        "maxTraceQ": args.max_trace_q,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "rows": rows,
        "interpretation": (
            "The closed trace condition alone is a nontrivial ODE constraint; "
            "local closed-trace jets can carry nonzero source rows.  The "
            "observability proof must use range inclusion/Volterra-Sturm "
            "structure, not bare ODE uniqueness."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
