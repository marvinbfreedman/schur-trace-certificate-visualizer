#!/usr/bin/env python3
r"""Finite noncollapse certificate for the projective response theorem.

For the fixed response map

    L f = (D_a^7 Lambda_a(f), D_a^8 Lambda_a(f))|_{a=.545},

the limiting noncollapse theorem is equivalent to saying that no active source
eigenvector lies in ``ker L``.  In a finite A-normalized high block this follows
from the stronger inequality

    max_{f in ker L, ||f||_A=1} <S f,f>  <  delta,

where ``delta`` is the source-active cutoff separating the two active
eigenvalues from the inactive tail.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import parse_orders, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, parse_ints  # noqa: E402
from lagrange_energy_control_certificate import split_spaces  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives, tower_row_on_polys  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from trace_to_source_kernel_refinement import active_indices  # noqa: E402


def response_matrix(args, base, e_derivs_point):
    out = mp.matrix(len(args.orders), base["scaledModes"].cols)
    for i, q in enumerate(args.orders):
        row = tower_row_on_polys(
            base["polys"],
            mp.mpf(args.point),
            e_derivs_point,
            q,
            args.jet_order,
        )
        on_high = row * base["scaledModes"]
        for j in range(on_high.cols):
            out[i, j] = on_high[0, j]
    return out


def eigvals(mat):
    if mat.rows == 0:
        return []
    return mp.eigsy(sym(mat), eigvals_only=True)


def max_eig(mat):
    vals = eigvals(mat)
    return vals[-1] if vals else mp.mpf("0")


def run_case(args, basis, e_derivs_at_s0):
    bargs = SimpleNamespace(**vars(args))
    bargs.basis = basis
    base = local_case(bargs, basis, e_derivs_at_s0)
    active_idx, active_floor = active_indices(
        base["sourceVals"], base["sourceTop"], mp.mpf(args.active_tol)
    )

    point_args = SimpleNamespace(**vars(args))
    point_args.basis = basis
    point_args.s0 = args.point
    _vals, e_derivs_point, _lam_derivs = exact_trace_derivatives(
        point_args, max(args.orders)
    )
    Lmat = response_matrix(args, base, e_derivs_point)
    response_gram = Lmat.T * Lmat
    N, _U, response_rank, response_nullity = split_spaces(Lmat, args.rank_tol)
    kernel_source_max = max_eig(N.T * base["sourceGram"] * N) if N.cols else mp.mpf("0")
    source_vals = list(base["sourceVals"])
    source_top = base["sourceTop"]
    active_min = min([source_vals[i] for i in active_idx], default=mp.mpf("0"))
    inactive_max = max(
        [source_vals[i] for i in range(len(source_vals)) if i not in active_idx],
        default=mp.mpf("0"),
    )
    cutoff_gap = active_floor - kernel_source_max
    min_active_gap = active_min - kernel_source_max if active_idx else mp.mpf("0")
    return {
        "basis": basis,
        "highModes": base["highModes"],
        "activeDim": len(active_idx),
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "sourceTop": f(source_top),
        "activeSourceEigenvalues": [f(source_vals[i]) for i in active_idx],
        "inactiveMaxSourceEigenvalue": f(inactive_max),
        "responseRank": response_rank,
        "responseNullity": response_nullity,
        "responseGramEigenvalues": [f(x) for x in eigvals(response_gram)],
        "kernelSourceMax": f(kernel_source_max),
        "kernelSourceFracOfTop": f(kernel_source_max / source_top if source_top else mp.mpf("0")),
        "cutoffGap": f(cutoff_gap),
        "cutoffGapFracOfTop": f(cutoff_gap / source_top if source_top else mp.mpf("0")),
        "minActiveGap": f(min_active_gap),
        "passesCutoffNoncollapse": bool(kernel_source_max < active_floor),
        "passesActiveMinNoncollapse": bool(active_idx and kernel_source_max < active_min),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--bases", default="18 20 22 24")
    parser.add_argument("--basis", type=int, default=24)
    parser.add_argument("--max-deriv", type=int, default=8)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-6")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=70)
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
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--json-out", default="noncollapse_kernel_source_gap.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    bases = parse_ints(args.bases)
    max_local_constraints = max(max(1, basis - args.local_constraint_offset) for basis in bases)
    max_q = max(args.max_trace_q, max_local_constraints - 1)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    print(
        f"Noncollapse kernel/source gap point={args.point} orders={args.orders} bases={bases}",
        flush=True,
    )
    print("  basis active_dim rank ker_source/top cutoff_gap/top pass", flush=True)
    rows = []
    for basis in bases:
        row = run_case(args, basis, e_derivs_at_s0)
        rows.append(row)
        print(
            f"  {basis:5d} {row['activeDim']:10d} {row['responseRank']:4d} "
            f"{fmt(mp.mpf(row['kernelSourceFracOfTop']), 8):>14} "
            f"{fmt(mp.mpf(row['cutoffGapFracOfTop']), 8):>14} "
            f"{row['passesCutoffNoncollapse']}",
            flush=True,
        )

    min_cutoff_gap = min(mp.mpf(row["cutoffGap"]) for row in rows)
    all_pass = all(row["passesCutoffNoncollapse"] for row in rows)
    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "point": f(mp.mpf(args.point)),
        "orders": args.orders,
        "bases": bases,
        "activeTol": f(mp.mpf(args.active_tol)),
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "allPassCutoffNoncollapse": all_pass,
        "minCutoffGap": f(min_cutoff_gap),
        "rows": rows,
        "interpretation": (
            "If kernelSourceMax is below the active spectral cutoff, no active "
            "source eigenvector can lie in ker L.  This is the finite operator "
            "form of the response noncollapse theorem."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
