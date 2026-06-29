#!/usr/bin/env python3
r"""Verify the trace Lagrange identity and size the adjoint source.

For

    P f = Lambda_s(f)=sum_{k=0}^8 a_k(s) f^(k)(s),
    a_k=e_k/k!,

the formal adjoint is

    P^* h = sum_{k=0}^8 (-1)^k D_s^k(a_k h).

The corrected Lagrange concomitant is

    B_P[h,f] =
      sum_{k=1}^8 sum_{j=0}^{k-1}
        (-1)^(k-1-j) D_s^(k-1-j)(a_k h) f^(j).

This script uses analytic derivatives of the confluent defect row e(s), then
checks

    D_s B_P[h,f] = h P f - f P^*h

on random jets.  It also reports P^*h_u and the endpoint-row norm, which are
the two residual objects that must be controlled by the commuted Sturm energy
after imposing P f=0.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_kb_confluent_mp import integrate  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature  # noqa: E402
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_concomitant_exact_derivatives import (  # noqa: E402
    center_taylor_mats,
    eigen_taylor,
)
from trace_concomitant_membership import (  # noqa: E402
    product_derivative,
    trace_green_concomitant_row,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def trace_row(e_derivs, jet_order: int):
    return [e_derivs[0][k] / mp.factorial(k) for k in range(jet_order)]


def trace_value(e_derivs, jets, jet_order: int):
    row = trace_row(e_derivs, jet_order)
    return mp.fsum(row[k] * jets[k] for k in range(jet_order))


def adjoint_source_value(e_derivs, h_derivs, jet_order: int):
    return mp.fsum(
        ((-1) ** k) * product_derivative(e_derivs, h_derivs, k, k)
        for k in range(jet_order)
    )


def concomitant_derivative_row(e_derivs, h_derivs, jet_order: int):
    """Return coefficients of D_s B_P[h,f] on f^(0)..f^(8)."""
    base = trace_green_concomitant_row(e_derivs, h_derivs, jet_order)
    out = [mp.mpf("0") for _ in range(jet_order)]

    # Differentiate the coefficient rows.
    for k in range(1, jet_order):
        for j in range(k):
            n = k - 1 - j
            out[j] += ((-1) ** n) * product_derivative(e_derivs, h_derivs, k, n + 1)

    # Differentiate f^(j).
    for j, coeff in enumerate(base):
        out[j + 1] += coeff
    return out


def row_norm(row):
    return mp.sqrt(mp.fsum(x * x for x in row))


def load_exact_trace(args):
    c = mp.pi
    s0 = mp.mpf(args.s0)
    big_order = args.jet_order + args.max_trace_q
    big, _segments = integrate(
        "kb",
        c,
        s0,
        big_order,
        mp.mpf(args.matrix_rmax),
        args.matrix_order,
    )
    mats = center_taylor_mats(big, args.jet_order, args.max_trace_q)
    vals, e_derivs, lam_derivs = eigen_taylor(mats, args.max_trace_q)
    return vals, e_derivs, lam_derivs


def random_jets(count: int, seed: int):
    rng = random.Random(seed)
    return [mp.mpf(rng.uniform(-1, 1)) for _ in range(count)]


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
    parser.add_argument("--random-tests", type=int, default=12)
    parser.add_argument("--seed", type=int, default=47)
    parser.add_argument("--json-out", default="trace_lagrange_adjoint_control.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    vals, e_derivs, lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    rows = []
    max_identity_defect = mp.mpf("0")
    max_relative_defect = mp.mpf("0")

    print(
        f"Trace Lagrange adjoint control s0={s0} max_trace_q={args.max_trace_q}",
        flush=True,
    )
    print(
        f"  lambda0={fmt(vals[0], 12)} gap={fmt(vals[1]-vals[0], 12)} "
        f"e8={fmt(e_derivs[0][-1], 12)}",
        flush=True,
    )
    print("  t       P*h          |B|          |P*h|/|B|    identity_defect", flush=True)

    for t_index, t in enumerate(t_values):
        h_derivs = source_derivatives(
            c,
            s0,
            t,
            args.jet_order,
            r_nodes,
            r_weights,
        )
        pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
        brow = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
        drow = concomitant_derivative_row(e_derivs, h_derivs, args.jet_order)
        tr = trace_row(e_derivs, args.jet_order)
        expected = [h_derivs[0] * tr[k] for k in range(args.jet_order)]
        expected[0] -= pstar
        row_defect = [drow[k] - expected[k] for k in range(args.jet_order)]
        row_defect_norm = row_norm(row_defect)

        random_max = mp.mpf("0")
        random_rel = mp.mpf("0")
        for test in range(args.random_tests):
            jets = random_jets(args.jet_order, args.seed + 1000 * t_index + test)
            lhs = mp.fsum(drow[k] * jets[k] for k in range(args.jet_order))
            rhs = h_derivs[0] * trace_value(e_derivs, jets, args.jet_order) - jets[0] * pstar
            defect = abs(lhs - rhs)
            scale = max(mp.mpf("1"), abs(lhs), abs(rhs))
            random_max = max(random_max, defect)
            random_rel = max(random_rel, defect / scale)

        max_identity_defect = max(max_identity_defect, row_defect_norm, random_max)
        max_relative_defect = max(max_relative_defect, random_rel)
        bnorm = row_norm(brow)
        ratio = abs(pstar) / bnorm if bnorm else mp.mpf("0")
        rows.append(
            {
                "t": f(t),
                "h0": f(h_derivs[0]),
                "pStarH": f(pstar),
                "boundaryRowNorm": f(bnorm),
                "pStarOverBoundaryNorm": f(ratio),
                "traceRowNorm": f(row_norm(tr)),
                "identityRowDefectNorm": f(row_defect_norm),
                "randomIdentityDefectMax": f(random_max),
                "randomIdentityRelativeDefectMax": f(random_rel),
            }
        )
        print(
            f"  {fmt(t, 6):>6} {fmt(pstar, 12):>12} {fmt(bnorm, 8):>12} "
            f"{fmt(ratio, 8):>12} {fmt(row_defect_norm, 8):>16}",
            flush=True,
        )

    data = {
        "s0": f(s0),
        "jetOrder": args.jet_order,
        "maxTraceQ": args.max_trace_q,
        "lambda0": f(vals[0]),
        "lambda1": f(vals[1]),
        "gap": f(vals[1] - vals[0]),
        "e8": f(e_derivs[0][-1]),
        "lambdaDerivatives": [f(x) for x in lam_derivs],
        "maxIdentityDefect": f(max_identity_defect),
        "maxIdentityRelativeDefect": f(max_relative_defect),
        "rows": rows,
        "nextControlTarget": (
            "Bound the endpoint pairing B_P[h_u,f]|_a^b and the integral "
            "of f P^*h_u by the commuted Sturm energy on the closed trace space."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  max identity defect={fmt(max_identity_defect, 12)}", flush=True)
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
