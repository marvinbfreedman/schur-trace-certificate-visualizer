#!/usr/bin/env python3
r"""Source-side full-theta tail certificate.

The previous full-theta check measured

    K_red(Phi_{<=N}) - K_red(tilde Phi_3)

as a reduced Volterra quadratic-form perturbation.  This script builds the
corresponding source-row object directly.

For a theta profile Psi, define the source row from

    h_u(s) = K_red^Psi(s,u).

The derivatives h_u^(k)(s0) are obtained by differentiating the reduced
Volterra integrand in the first variable using formal Taylor series.  Those
derivatives are then passed through the same Lagrange concomitant used in the
source/Riesz theorem:

    E_u^Psi f = (B_P[h_u,f](s0), P^*h_u(s0) f(s0)).

Finally

    S_Psi = int E_u^Psi* E_u^Psi du

is approximated by the same weighted source quadrature and A-normalized high
block used in ``source_side_riesz_rank_theorem.py``.  The literal finite tail is

    S_tail = S_{Phi_{<=N}} - S_{tilde Phi_3}.

The last section records an analytic exponential envelope for the omitted
theta modes n>N on the source window.  That envelope is intentionally stated as
a source-tail input bound; a fully formal proof still has to propagate it
through the Lagrange rows with interval arithmetic.
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

from determinant_gap_bound_diagnostic import eigvals_gram, parse_orders, submatrix_rows, sym  # noqa: E402
from global_trace_active_gap_scan import local_case  # noqa: E402
from global_trace_observability_gap import f, fmt, source_nodes  # noqa: E402
from lagrange_energy_control_certificate import boundary_functional_row, eval_functional_row  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from mp_partial_k_corr import weights_for  # noqa: E402
from quotient_factorization_mp import columns, endpoint_b_quadrature  # noqa: E402
from source_side_quadrature_refinement import weight_source_rows, weights_for_nodes  # noqa: E402
from source_side_riesz_rank_theorem import op_norm_rect, singular_values  # noqa: E402
from trace_active_derivative_rank import derivative_matrix  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value, trace_green_concomitant_row  # noqa: E402


def zero(n):
    return [mp.mpf("0") for _ in range(n)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scale(a, k):
    return [k * x for x in a]


def mul(a, b):
    n = min(len(a), len(b))
    out = zero(n)
    for i in range(n):
        out[i] = mp.fsum(a[j] * b[i - j] for j in range(i + 1))
    return out


def div(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = (
            a[i] - mp.fsum(b[j] * out[i - j] for j in range(1, i + 1))
        ) / b[0]
    return out


def exp_series(a):
    n = len(a)
    out = zero(n)
    out[0] = mp.e**a[0]
    for i in range(1, n):
        out[i] = mp.fsum(k * a[k] * out[i - k] for k in range(1, i + 1)) / i
    return out


def cosh_series(a):
    ep = exp_series(a)
    em = exp_series(scale(a, -1))
    return scale(add(ep, em), mp.mpf("0.5"))


def pieces_from_mode_weights(mode_weights):
    pieces = []
    for mode, weight in sorted(mode_weights.items()):
        n2 = mode * mode
        c = mp.pi * n2
        pieces.append((weight * 4 * c * c, mp.mpf("2.25"), c))
        pieces.append((weight * -6 * c, mp.mpf("1.25"), c))
    return pieces


def full_weights(nmax):
    return {n: mp.mpf("1") for n in range(1, nmax + 1)}


def psi_series(v0, pieces, n):
    ev0 = mp.e**v0
    e_delta = [mp.mpf(1) / mp.factorial(k) for k in range(n)]
    out = zero(n)
    for coeff, beta, c in pieces:
        exponent = scale(e_delta, -c * ev0)
        exponent[0] += beta * v0
        if n > 1:
            exponent[1] += beta
        out = add(out, scale(exp_series(exponent), coeff))
    return out


def psi_value(v, pieces):
    ev = mp.e**v
    return mp.fsum(coeff * mp.e ** (beta * v - c * ev) for coeff, beta, c in pieces)


def weight_series(center, omega, n):
    y = zero(n)
    y[0] = center
    if n > 1:
        y[1] = mp.mpf("0.5")
    return mul(y, cosh_series(scale(y, omega)))


def h_derivatives_for_pieces(args, pieces, u, r_nodes, r_weights):
    """Return h_u^(0)..h_u^(jet_order-1)(s0) for h_u(s)=K_red^Psi(s,u)."""
    n = args.jet_order
    s0 = mp.mpf(args.s0)
    omega = mp.mpf(args.omega)
    denom_s = psi_series(s0, pieces, n)
    denom_u = psi_value(u, pieces)
    ratio_u_cache = []
    for r in r_nodes:
        ratio_u_cache.append(psi_value(u + r, pieces) / denom_u)

    coeffs = zero(n)
    for r, weight, ratio_u in zip(r_nodes, r_weights, ratio_u_cache):
        ratio_s = div(psi_series(s0 + r, pieces, n), denom_s)
        w = weight_series(r + (s0 + u) / 2, omega, n)
        term = scale(mul(w, ratio_s), ratio_u)
        for k in range(n):
            coeffs[k] += weight * term[k]
    return [mp.factorial(k) * coeffs[k] for k in range(n)]


def source_rows_for_pieces(args, polys, e_derivs, pieces, r_nodes, r_weights, u):
    h_derivs = h_derivatives_for_pieces(args, pieces, u, r_nodes, r_weights)
    pstar = adjoint_source_value(e_derivs, h_derivs, args.jet_order)
    brow = trace_green_concomitant_row(e_derivs, h_derivs, args.jet_order)
    brow_poly = boundary_functional_row(polys, mp.mpf(args.s0), brow)
    eval_poly = eval_functional_row(polys, mp.mpf(args.s0), pstar)
    out = mp.matrix(2, len(polys))
    for j in range(len(polys)):
        out[0, j] = brow_poly[0, j]
        out[1, j] = eval_poly[0, j]
    return out


def weighted_source_matrix_for_pieces(args, base, e_derivs, pieces):
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.source_rmax), args.source_order)
    nodes = source_nodes(args)
    rows = mp.matrix(2 * len(nodes), len(base["polys"]))
    for block, u in enumerate(nodes):
        row = source_rows_for_pieces(
            args,
            base["polys"],
            e_derivs,
            pieces,
            r_nodes,
            r_weights,
            u,
        )
        for i in range(2):
            for j in range(len(base["polys"])):
                rows[2 * block + i, j] = row[i, j]
    weights = weights_for_nodes(nodes)
    weighted_rows = weight_source_rows(rows, weights)
    source_on_scaled = weighted_rows * base["scaledModes"]
    return nodes, weights, source_on_scaled, sym(source_on_scaled.T * source_on_scaled)


def op_norm_sym(mat):
    if mat.rows == 0:
        return mp.mpf("0")
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def eigensystem_source_gram(source_gram):
    vals, vecs = mp.eigsy(sym(source_gram), eigvals_only=False)
    if len(vals) < 3:
        raise ValueError("source Gram needs at least three eigenvalues")
    active_idx = [len(vals) - 2, len(vals) - 1]
    active = columns(vecs, active_idx)
    gap = vals[-2] - vals[-3]
    return vals, vecs, active_idx, active, gap


def finite_riesz_margin(args, base, active_basis):
    point_args = SimpleNamespace(**vars(args))
    point_args.s0 = args.point
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
    return lmat, finite_svals, rank_margin, response_norm


def derivative_tail_envelope(start_mode, vmin, max_order, stop_mode):
    """Conservative derivative envelope for omitted theta modes on v>=vmin."""
    ymin = mp.e**vmin
    per_order = [mp.mpf("0") for _ in range(max_order + 1)]
    per_mode = []
    for n in range(start_mode, stop_mode + 1):
        n2 = n * n
        c = mp.pi * n2
        mode_total = [mp.mpf("0") for _ in range(max_order + 1)]
        for coeff, beta in ((4 * c * c, mp.mpf("2.25")), (-6 * c, mp.mpf("1.25"))):
            base = abs(coeff) * ymin**beta * mp.e ** (-c * ymin)
            # For n>=9 and v>=0.08, c*exp(v) is much larger than beta+k.
            # The displayed bound uses |D^k exp(beta v-c exp v)| <=
            # exp(beta v-c exp v) (1+beta+c exp v)^k at the left endpoint.
            factor = 1 + beta + c * ymin
            for k in range(max_order + 1):
                mode_total[k] += base * factor**k
        for k in range(max_order + 1):
            per_order[k] += mode_total[k]
        if n in (start_mode, start_mode + 1, start_mode + 2, stop_mode):
            per_mode.append({"mode": n, "maxDerivativeBounds": [f(x) for x in mode_total]})
    return per_order, per_mode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=65)
    parser.add_argument("--full-nmax", type=int, default=8)
    parser.add_argument("--tail-start", type=int, default=9)
    parser.add_argument("--tail-stop", type=int, default=80)
    parser.add_argument("--source-order", type=int, default=28)
    parser.add_argument("--source-rmax", default="10")
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
    parser.add_argument("--dps", type=int, default=60)
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
    parser.add_argument("--json-out", default="full_theta_source_tail_certificate.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps
    if args.full_nmax < 4:
        raise SystemExit("--full-nmax must be at least 4")

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print(f"Preparing trace derivatives basis={args.basis} max_q={max_q}", flush=True)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, max_q)
    base_args = SimpleNamespace(**vars(args))
    base_args.source_grid = 9
    print("Preparing A-normalized high block", flush=True)
    base = local_case(base_args, args.basis, e_derivs_at_s0)

    tilde_pieces = pieces_from_mode_weights(weights_for("tilde3"))
    full_pieces = pieces_from_mode_weights(full_weights(args.full_nmax))

    print(f"Building S_tilde source grid={args.source_grid}", flush=True)
    nodes, weights, e_tilde, s_tilde = weighted_source_matrix_for_pieces(
        args, base, e_derivs_at_s0, tilde_pieces
    )
    print(f"Building S_full_{args.full_nmax}", flush=True)
    _nodes, _weights, e_full, s_full = weighted_source_matrix_for_pieces(
        args, base, e_derivs_at_s0, full_pieces
    )
    e_tail = e_full - e_tilde
    s_tail = sym(s_full - s_tilde)

    svals, _svecs, active_idx, active_basis, gap = eigensystem_source_gram(s_tilde)
    active_tail = sym(active_basis.T * s_tail * active_basis)
    active_core = sym(active_basis.T * s_tilde * active_basis)
    tail_norm = op_norm_sym(s_tail)
    core_norm = op_norm_sym(s_tilde)
    active_tail_norm = op_norm_sym(active_tail)
    active_core_norm = op_norm_sym(active_core)
    e_tail_norm = op_norm_rect(e_tail)
    e_core_norm = op_norm_rect(e_tilde)

    lmat, finite_svals, rank_margin, response_norm = finite_riesz_margin(args, base, active_basis)
    alpha_tail = 4 * tail_norm / gap if gap else mp.inf
    lower_after_tail = rank_margin * (1 - alpha_tail) - response_norm * alpha_tail

    vmin = min(mp.mpf(args.source_min), mp.mpf(args.s0))
    derivative_bounds, per_mode = derivative_tail_envelope(
        args.tail_start,
        vmin,
        args.jet_order - 1,
        args.tail_stop,
    )
    tail0 = derivative_bounds[0]
    psi_floor = min(abs(psi_value(x, tilde_pieces)) for x in [mp.mpf(args.source_min), mp.mpf(args.s0), mp.mpf(args.source_max)])
    relative_tail0 = tail0 / psi_floor if psi_floor else mp.inf

    print("Full-theta literal source-tail certificate", flush=True)
    print(f"  ||S_tail||/||S_tilde||={fmt(tail_norm/core_norm if core_norm else mp.inf, 12)}", flush=True)
    print(
        f"  active ||S_tail||/||S_tilde||={fmt(active_tail_norm/active_core_norm if active_core_norm else mp.inf, 12)}",
        flush=True,
    )
    print(f"  E-row tail/core={fmt(e_tail_norm/e_core_norm if e_core_norm else mp.inf, 12)}", flush=True)
    print(f"  gap={fmt(gap, 12)} alpha_tail={fmt(alpha_tail, 12)}", flush=True)
    print(f"  lower_after_tail={fmt(lower_after_tail, 12)}", flush=True)
    print(f"  n>={args.tail_start} derivative-0 envelope={fmt(tail0, 12)} relative={fmt(relative_tail0, 12)}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "sourceGrid": args.source_grid,
        "fullNmax": args.full_nmax,
        "sourceOrder": args.source_order,
        "sourceRmax": f(mp.mpf(args.source_rmax)),
        "sourceNodesFirstLast": [f(nodes[0]), f(nodes[-1])],
        "sourceWeightsFirstLast": [f(weights[0]), f(weights[-1])],
        "trace": {
            "lambda0AtS0": f(vals[0]),
            "lambda1AtS0": f(vals[1]),
            "gapAtS0": f(vals[1] - vals[0]),
            "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        },
        "sourceTailObject": {
            "definition": "S_tail = E_{Phi<=N}^* E_{Phi<=N} - E_tilde^* E_tilde",
            "coreOperatorNorm": f(core_norm),
            "tailOperatorNorm": f(tail_norm),
            "tailRelativeToCore": f(tail_norm / core_norm if core_norm else mp.inf),
            "rowTailOperatorNorm": f(e_tail_norm),
            "rowCoreOperatorNorm": f(e_core_norm),
            "rowTailRelativeToCore": f(e_tail_norm / e_core_norm if e_core_norm else mp.inf),
        },
        "activePlane": {
            "activeEigenvalues": [f(svals[i]) for i in active_idx],
            "complementTopEigenvalue": f(svals[-3]),
            "spectralGapToComplement": f(gap),
            "coreMatrix": [[f(active_core[i, j]) for j in range(active_core.cols)] for i in range(active_core.rows)],
            "tailMatrix": [[f(active_tail[i, j]) for j in range(active_tail.cols)] for i in range(active_tail.rows)],
            "tailOperatorNorm": f(active_tail_norm),
            "coreOperatorNorm": f(active_core_norm),
            "tailRelativeToCore": f(active_tail_norm / active_core_norm if active_core_norm else mp.inf),
        },
        "rieszStabilityUnderFiniteTail": {
            "tailProjectorAlpha": f(alpha_tail),
            "tailGapConditionPasses": bool(tail_norm < gap / 4),
            "finiteActiveResponseSingularValues": [f(x) for x in finite_svals],
            "finiteRankMargin": f(rank_margin),
            "responseOperatorNorm": f(response_norm),
            "lowerBoundAfterFiniteTail": f(lower_after_tail),
            "rankPersistsForPhiTruncation": bool(tail_norm < gap / 4 and lower_after_tail > 0),
        },
        "omittedThetaTailEnvelope": {
            "startMode": args.tail_start,
            "stopMode": args.tail_stop,
            "vMin": f(vmin),
            "yMin": f(mp.e**vmin),
            "psiFloorSampledOnSourceWindow": f(psi_floor),
            "maxDerivativeBounds": [f(x) for x in derivative_bounds],
            "relativeZerothDerivativeBound": f(relative_tail0),
            "sampledModes": per_mode,
            "interpretation": (
                "Conservative envelope for derivatives of the omitted theta profile on v>=v_min. "
                "This is the analytic input for the n>=tail_start source-row tail bound; the "
                "remaining formal step is interval propagation through the normalized source-row map."
            ),
        },
        "interpretation": (
            "Literal finite source-tail object for h_u(s)=K_red^Psi(s,u), built in the same "
            "A-normalized high-block coordinates as the source/Riesz theorem.  The finite "
            "Phi<=N tail is measured on the source Gram S=E^*E and through the Riesz gap/margin. "
            "The omitted n>=tail_start part is bounded at the theta-profile derivative level."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
