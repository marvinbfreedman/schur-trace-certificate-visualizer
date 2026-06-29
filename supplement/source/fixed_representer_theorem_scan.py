#!/usr/bin/env python3
r"""Refinement scan for the fixed-representer Hardy/Green theorem.

The closed-trace source-window theorem has been reduced to fixed Green
representers and scalar coefficient envelopes:

    adjoint eval:  p(u) k_0^{hi}
    boundary:      sum_{j=0}^7 b_j(u) k_j^{hi}.

The source coefficients p,p',p'',b,b',b'' are independent of the Galerkin
basis.  The basis-dependent part is only the fixed jet Gram matrix

    G_ij = <k_i^{hi}, k_j^{hi}>_A.

This script computes the scalar envelopes once and then scans nearby Galerkin
sections to check stability of the fixed-representer constants.
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

from adjoint_eval_representer_certificate import (  # noqa: E402
    energy_representer,
    source_derivatives_du_order,
)
from boundary_row_representer_certificate import (  # noqa: E402
    gram_norm,
    jet_eval_row,
    matrix_to_lists,
)
from lagrange_energy_control_certificate import make_qargs, split_spaces  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_b_quadrature,
    gram_matrix,
    trace_matrix,
)
from trace_lagrange_adjoint_control import (  # noqa: E402
    adjoint_source_value,
    load_exact_trace,
    trace_green_concomitant_row,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def constraints_for(args, basis: int) -> int:
    if args.constraint_rule == "fixed":
        return args.constraints
    if args.constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(args.constraint_ratio * basis)))
    if args.constraint_rule == "ratio-round":
        return max(1, min(basis - 1, round(args.constraint_ratio * basis)))
    if args.constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.constraint_offset))
    raise ValueError(args.constraint_rule)


def child_args(args, basis: int) -> SimpleNamespace:
    out = SimpleNamespace(**vars(args))
    out.basis = basis
    out.constraints = constraints_for(args, basis)
    return out


def coefficient_envelope(args):
    vals, e_derivs, _lam_derivs = load_exact_trace(args)
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)
    source_min = mp.mpf(args.source_min)
    source_max = mp.mpf(args.source_max)
    step = (source_max - source_min) / (args.source_grid - 1)
    radius = step / 2
    s0 = mp.mpf(args.s0)
    c = mp.pi

    rows = []
    max_abs_p = mp.mpf("0")
    max_abs_dp = mp.mpf("0")
    max_abs_d2p = mp.mpf("0")
    for i in range(args.source_grid):
        u = source_min + i * step
        h0 = source_derivatives_du_order(c, s0, u, args.jet_order, 0, r_nodes, r_weights)
        h1 = source_derivatives_du_order(c, s0, u, args.jet_order, 1, r_nodes, r_weights)
        h2 = source_derivatives_du_order(c, s0, u, args.jet_order, 2, r_nodes, r_weights)
        p = adjoint_source_value(e_derivs, h0, args.jet_order)
        dp = adjoint_source_value(e_derivs, h1, args.jet_order)
        d2p = adjoint_source_value(e_derivs, h2, args.jet_order)
        b = trace_green_concomitant_row(e_derivs, h0, args.jet_order)
        db = trace_green_concomitant_row(e_derivs, h1, args.jet_order)
        d2b = trace_green_concomitant_row(e_derivs, h2, args.jet_order)
        max_abs_p = max(max_abs_p, abs(p))
        max_abs_dp = max(max_abs_dp, abs(dp))
        max_abs_d2p = max(max_abs_d2p, abs(d2p))
        rows.append(
            {
                "u": f(u),
                "p": f(p),
                "dp": f(dp),
                "d2p": f(d2p),
                "b": [f(x) for x in b],
                "db": [f(x) for x in db],
                "d2b": [f(x) for x in d2b],
            }
        )
    return {
        "traceLambda0": f(vals[0]),
        "traceLambda1": f(vals[1]),
        "sourceMin": f(source_min),
        "sourceMax": f(source_max),
        "sourceGrid": args.source_grid,
        "mesh": f(step),
        "meshRadius": f(radius),
        "maxAbsP": f(max_abs_p),
        "maxAbsDP": f(max_abs_dp),
        "maxAbsD2P": f(max_abs_d2p),
        "coverAbsP": f(max_abs_p + radius * max_abs_dp),
        "coverAbsDP": f(max_abs_dp + radius * max_abs_d2p),
        "rows": rows,
    }


def fixed_jet_gram(args, coeff_env, basis: int):
    bargs = child_args(args, basis)
    qargs = make_qargs(bargs)
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

    order = args.jet_order - 1
    reps = []
    jet_records = []
    for deriv in range(order):
        row_n = jet_eval_row(polys, mp.mpf(args.s0), deriv) * N
        rep_n, norm2, rel = energy_representer(row_n, A, high_a_modes, high_avals)
        reps.append(rep_n)
        jet_records.append(
            {
                "deriv": deriv,
                "norm2": f(norm2),
                "norm": f(mp.sqrt(max(mp.mpf("0"), norm2))),
                "rangeRelativeDefect": f(rel),
            }
        )

    G = mp.matrix(order)
    for i in range(order):
        for j in range(order):
            G[i, j] = (reps[i].T * A * reps[j])[0, 0]

    radius = mp.mpf(coeff_env["meshRadius"])
    cover_abs_p = mp.mpf(coeff_env["coverAbsP"])
    cover_abs_dp = mp.mpf(coeff_env["coverAbsDP"])
    eval_norm2 = G[0, 0]
    cover_eval = cover_abs_p * cover_abs_p * eval_norm2
    cover_deval = cover_abs_dp * cover_abs_dp * eval_norm2

    max_b_norm = mp.mpf("0")
    max_db_norm = mp.mpf("0")
    max_d2b_norm = mp.mpf("0")
    for row in coeff_env["rows"]:
        b = [mp.mpf(x) for x in row["b"]]
        db = [mp.mpf(x) for x in row["db"]]
        d2b = [mp.mpf(x) for x in row["d2b"]]
        max_b_norm = max(max_b_norm, gram_norm(G, b))
        max_db_norm = max(max_db_norm, gram_norm(G, db))
        max_d2b_norm = max(max_d2b_norm, gram_norm(G, d2b))

    cover_b = max_b_norm + radius * max_db_norm
    cover_db = max_db_norm + radius * max_d2b_norm
    return {
        "basis": basis,
        "constraints": bargs.constraints,
        "rank": rank,
        "nullity": nullity,
        "positiveModes": len(avals),
        "cutoff": cutoff,
        "highModes": len(high_avals),
        "lambdaMinHigh": f(high_avals[0]) if high_avals else None,
        "lambdaMaxHigh": f(high_avals[-1]) if high_avals else None,
        "evalRepresenterNorm2": f(eval_norm2),
        "evalCoverNorm2": f(cover_eval),
        "evalDerivativeCoverNorm2": f(cover_deval),
        "maxBNormGrid": f(max_b_norm),
        "maxDBNormGrid": f(max_db_norm),
        "maxD2BNormGrid": f(max_d2b_norm),
        "boundaryCoverNorm2": f(cover_b * cover_b),
        "boundaryDerivativeCoverNorm2": f(cover_db * cover_db),
        "maxJetRangeRelativeDefect": f(max(mp.mpf(r["rangeRelativeDefect"]) for r in jet_records)),
        "jetRepresenters": jet_records,
        "jetGram": matrix_to_lists(G),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=17)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20 22")
    parser.add_argument("--basis", type=int, default=20)
    parser.add_argument("--constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--constraint-offset", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
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
    parser.add_argument("--json-out", default="fixed_representer_theorem_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    print(
        f"Fixed-representer theorem scan bases={bases} source_grid={args.source_grid}",
        flush=True,
    )
    coeff_env = coefficient_envelope(args)
    print(
        f"  scalar covers |p|={fmt(mp.mpf(coeff_env['coverAbsP']), 10)} "
        f"|p'|={fmt(mp.mpf(coeff_env['coverAbsDP']), 10)}",
        flush=True,
    )
    print("  basis cons pos eval_norm2 eval_cover bdry_cover dB_cover", flush=True)
    results = []
    for basis in bases:
        row = fixed_jet_gram(args, coeff_env, basis)
        results.append(row)
        print(
            f"  {basis:5d} {row['constraints']:4d} {row['positiveModes']:3d} "
            f"{fmt(mp.mpf(row['evalRepresenterNorm2']), 8):>11} "
            f"{fmt(mp.mpf(row['evalCoverNorm2']), 8):>11} "
            f"{fmt(mp.mpf(row['boundaryCoverNorm2']), 8):>11} "
            f"{fmt(mp.mpf(row['boundaryDerivativeCoverNorm2']), 8):>11}",
            flush=True,
        )

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "bases": bases,
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "cutoff": args.cutoff,
        "coefficientEnvelope": coeff_env,
        "rows": results,
        "theoremTemplate": (
            "Fixed-representer Hardy/Green theorem. If the closed-trace high "
            "RKHS has bounded jet representers k_j^{hi}, and the source "
            "coefficient rows p,p',b,b' are bounded on the compact source "
            "window, then the source-window Hardy/Green inequality follows."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
