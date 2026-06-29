#!/usr/bin/env python3
r"""Riesz-projector rank theorem from the rigorous source quadrature bound.

This closes the source-side noncollapse step in the finite-core program.
The continuum source Gram is

    S = E^* E = int E(u)^* E(u) du

on the fixed A-normalized high block.  The complex-ball Bernstein certificate
gives an absolute operator bound

    ||S - S_N|| <= eps

for the weighted trapezoid source Gram ``S_N``.  If the top two eigenvalues of
``S_N`` are separated from the rest by ``g`` and ``eps < g/4``, the Riesz
projector ``P`` for the continuum top-two cluster obeys the conservative bound

    ||P-P_N|| <= alpha := 4 eps/g.

Let ``m = sigma_min(L P_N)`` on the finite active two-plane and let ``ell=||L||``.
For every unit vector ``v`` in Range(P),

    ||L v|| >= m(1-alpha) - ell alpha.

Thus ``L`` has rank two on the continuum active right eigenspace.  Since the
positive eigenvectors of ``T=E E^*`` map to those of ``S=E^*E`` by ``E^*``,
this proves

    rank(L E^* | U_delta) = 2.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from determinant_gap_bound_diagnostic import eigvals_gram, parse_orders, submatrix_rows, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, stack_source_rows  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import columns  # noqa: E402
from source_side_quadrature_refinement import weight_source_rows, weights_for_nodes  # noqa: E402
from trace_active_derivative_rank import derivative_matrix  # noqa: E402


def load_quadrature_bound(path: str):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return data, mp.mpf(data["sourceGramErrorBound"])


def op_norm_rect(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = eigvals_gram(mat)
    return mp.sqrt(max([mp.mpf("0"), *vals]))


def singular_values(mat):
    vals = eigvals_gram(mat)
    return [mp.sqrt(max(mp.mpf("0"), val)) for val in vals]


def finite_weighted_source_gram(args, base, e_derivs, grid):
    gargs = SimpleNamespace(**vars(args))
    gargs.source_grid = grid
    nodes, source_rows = stack_source_rows(gargs, base["polys"], e_derivs)
    weights = weights_for_nodes(nodes)
    weighted_rows = weight_source_rows(source_rows, weights)
    source_on_scaled = weighted_rows * base["scaledModes"]
    return nodes, weights, sym(source_on_scaled.T * source_on_scaled), source_on_scaled


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quadrature-json", default="source_quadrature_complex_ball_bernstein_rho2.json")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=129)
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
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
    parser.add_argument("--json-out", default="source_side_riesz_rank_theorem.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps

    quad_data, eps = load_quadrature_bound(args.quadrature_json)
    if int(quad_data["sourceGrid"]) != args.source_grid:
        raise SystemExit(
            f"quadrature JSON grid {quad_data['sourceGrid']} does not match --source-grid {args.source_grid}"
        )

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace data max_q={max_q} basis={args.basis}", flush=True)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)

    # The source grid does not affect the A-normalized high block.  Use a small
    # temporary grid inside local_case to avoid building an unused dense source
    # sample before constructing the actual weighted S_N below.
    base_args = SimpleNamespace(**vars(args))
    base_args.source_grid = 9
    print("Preparing local high block", flush=True)
    base = local_case(base_args, args.basis, e_derivs_at_s0)

    print(f"Building weighted S_N grid={args.source_grid}", flush=True)
    nodes, weights, source_gram, source_on_scaled = finite_weighted_source_gram(
        args,
        base,
        e_derivs_at_s0,
        args.source_grid,
    )
    svals, svecs = mp.eigsy(source_gram, eigvals_only=False)
    if len(svals) < 3:
        raise SystemExit("source Gram has fewer than three eigenvalues")
    active_idx = [len(svals) - 2, len(svals) - 1]
    complement_top_idx = len(svals) - 3
    active_basis = columns(svecs, active_idx)
    source_top = svals[-1]
    active_min = svals[-2]
    complement_top = svals[complement_top_idx]
    gap = active_min - complement_top

    point_args = SimpleNamespace(**vars(args))
    point_args.s0 = args.point
    _vals_p, _e_derivs_point, _lam_derivs_p = exact_trace_derivatives(point_args, max(args.orders))
    full_response = derivative_matrix(
        point_args,
        base,
        base["scaledModes"],
        mp.mpf(args.point),
        max(args.orders),
    )
    lmat = submatrix_rows(full_response, args.orders)
    finite_active_response = lmat * active_basis
    finite_svals = singular_values(finite_active_response)
    rank_margin = min(finite_svals)
    response_norm = op_norm_rect(lmat)
    alpha = 4 * eps / gap if gap else mp.inf
    projector_gap_pass = bool(eps < gap / 4)
    lower_bound = rank_margin * (1 - alpha) - response_norm * alpha
    rank_pass = bool(projector_gap_pass and lower_bound > 0)

    active_floor = mp.mpf(args.active_tol) * max(mp.mpf("1"), source_top)
    active_dim_by_tol = sum(1 for val in svals if val > active_floor)
    source_norm = op_norm_rect(source_on_scaled)
    source_gram_norm = source_norm**2
    eps_relative_top = eps / source_top if source_top else mp.inf
    gap_relative_top = gap / source_top if source_top else mp.inf
    rank_margin_relative = rank_margin / response_norm if response_norm else mp.inf
    alpha_allowed = rank_margin / (rank_margin + response_norm) if rank_margin + response_norm else mp.mpf("0")

    print("Riesz rank theorem certificate", flush=True)
    print(f"  eps={fmt(eps, 12)}", flush=True)
    print(f"  source_top={fmt(source_top, 12)} active_min={fmt(active_min, 12)}", flush=True)
    print(f"  gap={fmt(gap, 12)} eps/gap={fmt(eps/gap, 12)} alpha={fmt(alpha, 12)}", flush=True)
    print(f"  margin={fmt(rank_margin, 12)} ||L||={fmt(response_norm, 12)}", flush=True)
    print(f"  lower_bound={fmt(lower_bound, 12)} pass={rank_pass}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "s0": f(mp.mpf(args.s0)),
        "point": f(mp.mpf(args.point)),
        "orders": args.orders,
        "basis": args.basis,
        "sourceMin": f(mp.mpf(args.source_min)),
        "sourceMax": f(mp.mpf(args.source_max)),
        "sourceGrid": args.source_grid,
        "quadratureJson": args.quadrature_json,
        "quadratureMethod": quad_data.get("interpretation", ""),
        "sourceGramErrorBound": f(eps),
        "sourceGramErrorRelativeToTop": f(eps_relative_top),
        "sourceGramNormFromFiniteGrid": f(source_gram_norm),
        "sourceTop": f(source_top),
        "activeEigenvalues": [f(svals[i]) for i in active_idx],
        "complementTopEigenvalue": f(complement_top),
        "spectralGapToComplement": f(gap),
        "spectralGapToComplementFracOfTop": f(gap_relative_top),
        "activeTol": f(mp.mpf(args.active_tol)),
        "activeFloor": f(active_floor),
        "activeDimByTol": active_dim_by_tol,
        "projectorErrorAlpha": f(alpha),
        "projectorErrorAlphaAllowed": f(alpha_allowed),
        "projectorGapConditionPasses": projector_gap_pass,
        "finiteActiveResponseSingularValues": [f(x) for x in finite_svals],
        "finiteRankMargin": f(rank_margin),
        "responseOperatorNorm": f(response_norm),
        "relativeRankMargin": f(rank_margin_relative),
        "continuumLowerBoundForUnitActiveVector": f(lower_bound),
        "rankTheoremPasses": rank_pass,
        "sourceNodesFirstLast": [f(nodes[0]), f(nodes[-1])],
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "theorem": (
            "For S=E^*E, ||S-S_N||<=eps and eps<gap/4 imply "
            "||P-P_N||<=4eps/gap for the active top-two Riesz projectors. "
            "Since sigma_min(LP_N)=m and ||L||=ell, every unit v in Range(P) "
            "satisfies ||Lv||>=m(1-alpha)-ell alpha.  The displayed positive "
            "lower bound proves rank(L on Range(P))=2, equivalently "
            "rank(LE^*|U_delta)=2 for the positive source eigenspace."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
