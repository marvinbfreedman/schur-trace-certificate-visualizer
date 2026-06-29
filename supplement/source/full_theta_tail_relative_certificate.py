#!/usr/bin/env python3
r"""Finite reduced-Volterra certificate for the full theta tail.

The finite-core theorem is built for the zero-slope corrected core

    \tilde\Phi_3 = \phi_1 + \phi_2 + \alpha_3 \phi_3.

The full kernel uses

    \Phi = \sum_{n>=1} \phi_n,

so the perturbation is

    \Phi - \tilde\Phi_3 = (1-\alpha_3)\phi_3 + \sum_{n>=4}\phi_n.

This script measures the induced perturbation in the same reduced Volterra
Galerkin coordinates used by the quotient/source program:

    K_tail = K_red(\Phi_{<=N}) - K_red(\tilde\Phi_3).

It also projects that tail onto the active source plane from the recently
certified Riesz theorem.  This is not yet the literal analytic ``S_tail`` for
the continuum source map; it is the finite normalized operator check that tells
us whether the theta tail is perturbative in the coordinates that matter.
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
from global_trace_observability_gap import f, fmt  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from mp_partial_k_corr import mode_right_derivative, weights_for  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_ratios,
    laguerre_integral,
    legendre_quadrature,
    shifted_legendre_polys,
    trace_matrix,
)
from source_side_quadrature_refinement import weight_source_rows, weights_for_nodes  # noqa: E402
from source_side_riesz_rank_theorem import finite_weighted_source_gram, op_norm_rect  # noqa: E402


def pieces_from_mode_weights(weights):
    pieces = []
    for mode, weight in sorted(weights.items()):
        n2 = mode * mode
        c = mp.pi * n2
        pieces.append((weight * 4 * c * c, mp.mpf("2.25"), c))
        pieces.append((weight * -6 * c, mp.mpf("1.25"), c))
    return pieces


def full_weights(nmax):
    return {n: mp.mpf("1") for n in range(1, nmax + 1)}


def tail_weights_from_tilde(nmax):
    tilde = weights_for("tilde3")
    out = {}
    for n in range(1, nmax + 1):
        value = mp.mpf("1") - tilde.get(n, mp.mpf("0"))
        if value:
            out[n] = value
    return out


def gram_from_pieces(args, pieces, polys):
    pts, wts = legendre_quadrature(mp.mpf(args.L), args.quad)
    lag_nodes, lag_weights = mp.gauss_quadrature(args.laguerre, "laguerre")
    values = mp.matrix(args.quad, args.basis)
    for i, x in enumerate(pts):
        for j, poly in enumerate(polys):
            total = mp.mpf("0")
            for coeff in reversed(poly):
                total = total * x + coeff
            values[i, j] = total

    ratio_cache = [endpoint_ratios(s, pieces) for s in pts]
    omega = mp.mpf(args.omega)

    def kernel_from_cached(a_idx, b_idx):
        s = pts[a_idx]
        t = pts[b_idx]
        center = (s + t) / 2
        es = mp.e**s
        et = mp.e**t
        total = mp.mpf("0")
        for ratio_i, beta_i, c_i in ratio_cache[a_idx]:
            for ratio_j, beta_j, c_j in ratio_cache[b_idx]:
                alpha = c_i * es + c_j * et
                p = beta_i + beta_j - 1
                total += ratio_i * ratio_j * laguerre_integral(
                    alpha,
                    p,
                    center,
                    omega,
                    lag_nodes,
                    lag_weights,
                )
        return total

    out = mp.matrix(args.basis)
    for a in range(args.quad):
        for b in range(a, args.quad):
            kval = kernel_from_cached(a, b)
            weight = wts[a] * wts[b] * kval
            if a == b:
                for i in range(args.basis):
                    vi = values[a, i]
                    for j in range(args.basis):
                        out[i, j] += weight * vi * values[b, j]
            else:
                for i in range(args.basis):
                    vai = values[a, i]
                    vbi = values[b, i]
                    for j in range(args.basis):
                        out[i, j] += weight * (
                            vai * values[b, j] + vbi * values[a, j]
                        )
    return sym(out)


def eig_op_norm_symmetric(mat):
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    return max([abs(v) for v in vals] + [mp.mpf("0")])


def positive_min_eig(mat, tol):
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    positives = [v for v in vals if v > tol]
    return positives[0] if positives else mp.mpf("0")


def split_trace(polys, args):
    centers, R = trace_matrix(polys, args)
    gram = R.T * R
    vals, vecs = mp.eigsy(sym(gram), eigvals_only=False)
    rmax = max([abs(v) for v in vals] + [mp.mpf("0")])
    tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    null_idx = [i for i, v in enumerate(vals) if v <= tol]
    range_idx = [i for i, v in enumerate(vals) if v > tol]
    return centers, R, columns(vecs, null_idx), columns(vecs, range_idx), tol


def project_active_plane(args, polys, e_derivs):
    # Reuse exactly the endpoint/source active plane used by the Riesz theorem.
    base_args = SimpleNamespace(**vars(args))
    base_args.source_grid = 9
    base = local_case(base_args, args.basis, e_derivs)
    nodes, weights, source_gram, _source_on_scaled = finite_weighted_source_gram(
        args,
        base,
        e_derivs,
        args.source_grid,
    )
    svals, svecs = mp.eigsy(sym(source_gram), eigvals_only=False)
    if len(svals) < 2:
        raise ValueError("source Gram has fewer than two eigenvalues")
    active_idx = [len(svals) - 2, len(svals) - 1]
    active_source = columns(svecs, active_idx)
    active_poly = base["scaledModes"] * active_source
    return {
        "nodes": nodes,
        "weights": weights,
        "sourceEigenvalues": [svals[i] for i in active_idx],
        "sourceTop": svals[-1],
        "sourceGap": svals[-2] - svals[-3] if len(svals) >= 3 else mp.mpf("0"),
        "activePoly": active_poly,
    }


def matrix_stats(label, mat, core=None):
    op = eig_op_norm_symmetric(mat)
    vals = mp.eigsy(sym(mat), eigvals_only=True)
    row = {
        f"{label}OperatorNorm": f(op),
        f"{label}MinEigenvalue": f(vals[0]),
        f"{label}MaxEigenvalue": f(vals[-1]),
    }
    if core is not None:
        core_op = eig_op_norm_symmetric(core)
        row[f"{label}RelativeToCoreOp"] = f(op / core_op if core_op else mp.inf)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=22)
    parser.add_argument("--quad", type=int, default=32)
    parser.add_argument("--laguerre", type=int, default=36)
    parser.add_argument("--full-nmax", type=int, default=12)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=129)
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--point", default="0.545")
    parser.add_argument("--orders", default="7 8")
    parser.add_argument("--max-deriv", type=int, default=8)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=12)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--active-tol", default="1e-6")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--riesz-json", default="source_side_riesz_rank_theorem.json")
    parser.add_argument("--json-out", default="full_theta_tail_relative_certificate.json")
    args = parser.parse_args()

    args.orders = parse_orders(args.orders)
    mp.mp.dps = args.dps

    if args.full_nmax < 4:
        raise SystemExit("--full-nmax must be at least 4")

    riesz = json.loads(Path(args.riesz_json).read_text(encoding="utf-8"))
    certified_lower = mp.mpf(riesz["continuumLowerBoundForUnitActiveVector"])
    finite_rank_margin = mp.mpf(riesz["finiteRankMargin"])
    projector_alpha = mp.mpf(riesz["projectorErrorAlpha"])

    alpha3 = weights_for("tilde3")[3]
    tilde_weights = weights_for("tilde3")
    full_w = full_weights(args.full_nmax)
    tail_w = tail_weights_from_tilde(args.full_nmax)
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))

    print(
        f"Full theta tail relative certificate basis={args.basis} "
        f"quad={args.quad} nmax={args.full_nmax}",
        flush=True,
    )
    print("  building K_tilde", flush=True)
    k_tilde = gram_from_pieces(args, pieces_from_mode_weights(tilde_weights), polys)
    print("  building K_full_truncated", flush=True)
    k_full = gram_from_pieces(args, pieces_from_mode_weights(full_w), polys)
    k_tail = sym(k_full - k_tilde)

    centers, R, N, U, trace_tol = split_trace(polys, args)
    k_tilde_ker = sym(N.T * k_tilde * N) if N.cols else mp.matrix(0)
    k_tail_ker = sym(N.T * k_tail * N) if N.cols else mp.matrix(0)
    k_tilde_range = sym(U.T * k_tilde * U) if U.cols else mp.matrix(0)
    k_tail_range = sym(U.T * k_tail * U) if U.cols else mp.matrix(0)

    max_q = max(args.max_trace_q, max(1, args.basis - args.local_constraint_offset) - 1)
    print("  building active source plane", flush=True)
    vals, e_derivs, lam_derivs = exact_trace_derivatives(args, max_q)
    active = project_active_plane(args, polys, e_derivs)
    V = active["activePoly"]
    active_core = sym(V.T * k_tilde * V)
    active_tail = sym(V.T * k_tail * V)
    active_full = sym(V.T * k_full * V)

    tail_active_norm = eig_op_norm_symmetric(active_tail)
    active_core_norm = eig_op_norm_symmetric(active_core)
    active_tail_relative = tail_active_norm / active_core_norm if active_core_norm else mp.inf

    global_tail_norm = eig_op_norm_symmetric(k_tail)
    global_core_norm = eig_op_norm_symmetric(k_tilde)
    ker_tail_norm = eig_op_norm_symmetric(k_tail_ker) if N.cols else mp.mpf("0")
    ker_core_min = positive_min_eig(k_tilde_ker, mp.mpf(args.psd_tol)) if N.cols else mp.mpf("0")
    trace_tail_norm = eig_op_norm_symmetric(k_tail_range) if U.cols else mp.mpf("0")

    # Unit warning: the Riesz rank lower bound is an L-response margin, while
    # K_tail is a reduced Volterra quadratic-form perturbation.  This ratio is
    # only a scale diagnostic, not the final analytic inequality.
    diagnostic_margin_ratio = tail_active_norm / certified_lower if certified_lower else mp.inf
    finite_margin_ratio = tail_active_norm / finite_rank_margin if finite_rank_margin else mp.inf
    perturbative_pass = bool(active_tail_relative < mp.mpf("1e-6"))

    print(f"  alpha3={fmt(alpha3, 18)}", flush=True)
    print(f"  global ||K_tail||/||K_core||={fmt(global_tail_norm / global_core_norm, 12)}", flush=True)
    print(f"  active ||tail||/||core||={fmt(active_tail_relative, 12)}", flush=True)
    print(f"  active tail / Riesz lower diagnostic={fmt(diagnostic_margin_ratio, 12)}", flush=True)

    data = {
        "omega": f(mp.mpf(args.omega)),
        "L": f(mp.mpf(args.L)),
        "basis": args.basis,
        "quad": args.quad,
        "laguerre": args.laguerre,
        "fullNmax": args.full_nmax,
        "sourceGrid": args.source_grid,
        "alpha3": f(alpha3),
        "tildeWeights": {str(k): f(v) for k, v in tilde_weights.items()},
        "tailWeights": {str(k): f(v) for k, v in tail_w.items()},
        "traceRank": U.cols,
        "traceNullity": N.cols,
        "traceRankTolerance": f(trace_tol),
        "traceCenters": [f(x) for x in centers],
        "global": {
            **matrix_stats("core", k_tilde),
            **matrix_stats("tail", k_tail, k_tilde),
            **matrix_stats("full", k_full),
        },
        "traceKernel": {
            "coreKerPositiveMinEigenvalue": f(ker_core_min),
            "tailKerOperatorNorm": f(ker_tail_norm),
            "tailKerRelativeToCorePositiveMin": f(ker_tail_norm / ker_core_min if ker_core_min else mp.inf),
            "tailTraceRangeOperatorNorm": f(trace_tail_norm),
        },
        "activeSourcePlane": {
            "sourceEigenvalues": [f(x) for x in active["sourceEigenvalues"]],
            "sourceTop": f(active["sourceTop"]),
            "sourceGap": f(active["sourceGap"]),
            "coreMatrix": [[f(active_core[i, j]) for j in range(active_core.cols)] for i in range(active_core.rows)],
            "tailMatrix": [[f(active_tail[i, j]) for j in range(active_tail.cols)] for i in range(active_tail.rows)],
            "fullMatrix": [[f(active_full[i, j]) for j in range(active_full.cols)] for i in range(active_full.rows)],
            "tailOperatorNorm": f(tail_active_norm),
            "coreOperatorNorm": f(active_core_norm),
            "tailRelativeToCore": f(active_tail_relative),
            "tailOverCertifiedRieszLowerDiagnostic": f(diagnostic_margin_ratio),
            "tailOverFiniteRankMarginDiagnostic": f(finite_margin_ratio),
            "projectorAlphaFromRieszCertificate": f(projector_alpha),
            "certifiedRieszLowerBound": f(certified_lower),
            "finiteRankMargin": f(finite_rank_margin),
            "perturbativeDiagnosticPasses": perturbative_pass,
        },
        "lambda0AtS0": f(vals[0]),
        "lambda1AtS0": f(vals[1]),
        "gapAtS0": f(vals[1] - vals[0]),
        "lambdaDerivativesAtS0": [f(x) for x in lam_derivs],
        "interpretation": (
            "Finite reduced-Volterra theta-tail diagnostic in the same "
            "polynomial/source-active coordinates as the Riesz theorem.  This "
            "does not yet replace the literal continuum S_tail derivation; it "
            "checks that K_red(Phi)-K_red(tildePhi_3) is small on the active "
            "source plane and trace-normalized Galerkin blocks."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
