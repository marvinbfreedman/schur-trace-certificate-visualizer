#!/usr/bin/env python3
r"""Finite construction attempt for z -> x_z using only mpmath.

We try to construct trace coordinates x_z in the completed Volterra trace
model so that the Volterra/KLM branch Grams match the Hardy branch Grams.

The finite method is conservative:

  1. Build the finite Green-minimizer trace image.
  2. Build Volterra plus/minus/cross Gram forms P, M, C on trace coordinates.
  3. Choose X=(x_z) so X^* P X matches the Hardy plus Gram exactly.
  4. Search the remaining orthonormal freedom for the best joint Gram match.

This can prove a finite candidate if the joint residual is tiny.  Otherwise it
identifies the missing ingredient: x_z cannot be chosen by plus-Gram matching
alone; it must come from the actual trace functionals applied to the Hardy
evaluation exponential.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import (  # noqa: E402
    default_z_nodes,
    hardy_feature_scalar,
    halfline_inner,
)
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    endpoint_ratios,
    frob_norm,
    poly_value,
    positive_part_inverse,
    shifted_legendre_polys,
    trace_matrix,
)
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


def ctranspose(a: mp.matrix) -> mp.matrix:
    out = mp.matrix(a.cols, a.rows)
    for i in range(a.rows):
        for j in range(a.cols):
            out[j, i] = mp.conj(a[i, j])
    return out


def herm(a: mp.matrix) -> mp.matrix:
    return (a + ctranspose(a)) / 2


def sym(a: mp.matrix) -> mp.matrix:
    return (a + a.T) / 2


def block2(a: mp.matrix, b: mp.matrix, c: mp.matrix, d: mp.matrix) -> mp.matrix:
    out = mp.matrix(a.rows + d.rows, a.cols + d.cols)
    for i in range(a.rows):
        for j in range(a.cols):
            out[i, j] = a[i, j]
    for i in range(b.rows):
        for j in range(b.cols):
            out[i, a.cols + j] = b[i, j]
    for i in range(c.rows):
        for j in range(c.cols):
            out[a.rows + i, j] = c[i, j]
    for i in range(d.rows):
        for j in range(d.cols):
            out[a.rows + i, a.cols + j] = d[i, j]
    return out


def rel_frob(a: mp.matrix, b: mp.matrix) -> mp.mpf:
    return frob_norm(a - b) / max(frob_norm(b), mp.mpf("1e-300"))


def local_args(args: argparse.Namespace):
    return SimpleNamespace(
        constraints=args.constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
    )


def legendre_quadrature(length: mp.mpf, order: int):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    pts = [mp.mpf("0.5") * length * (x + 1) for x in nodes]
    wts = [mp.mpf("0.5") * length * w for w in weights]
    return pts, wts


def a_value(s: mp.mpf, u: mp.mpf, pcs) -> mp.mpf:
    es = mp.e**s
    eu = mp.e**u
    total = mp.mpf("0")
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * mp.e ** (beta * u - c * es * (eu - 1))
    return total


def branch_feature_matrices_mp(args: argparse.Namespace, sigma: int):
    pcs = pieces_for(args.kind)
    s_pts, s_wts = legendre_quadrature(mp.mpf(args.L), args.s_order)
    u_pts, u_wts = legendre_quadrature(mp.mpf(args.u_max), args.u_order)
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    values = mp.matrix(len(s_pts), args.basis)
    for i, s in enumerate(s_pts):
        for j, poly in enumerate(polys):
            values[i, j] = poly_value(poly, s)
    m = mp.matrix(args.u_order, args.basis)
    n = mp.matrix(args.u_order, args.basis)
    for iu, u in enumerate(u_pts):
        root_u = mp.sqrt(u_wts[iu])
        for js, s in enumerate(s_pts):
            branch = mp.e ** (mp.mpf("0.5") * sigma * mp.mpf(args.omega) * (s + u))
            base = root_u * s_wts[js] * a_value(s, u, pcs) * branch
            for k in range(args.basis):
                row_value = base * values[js, k]
                m[iu, k] += row_value
                n[iu, k] += (s + u) * row_value
    return m, n, u_pts


def quotient_minimizer(K: mp.matrix, R: mp.matrix, args: argparse.Namespace):
    gram = R.T * R
    rvals, rvecs = mp.eigsy(sym(gram), eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)
    A = N.T * K * N
    B = N.T * K * U
    _vals, _vecs, _keep, _zero, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    return U - N * Aplus * B, U.cols, N.cols, U


def finite_trace_model(args: argparse.Namespace):
    sigmas = [0] if mp.mpf(args.omega) == 0 else [-1, 1]
    branch_weight = mp.mpf("1") if mp.mpf(args.omega) == 0 else mp.mpf("0.5")
    pieces = []
    K = mp.matrix(args.basis)
    for sigma in sigmas:
        m, n, u_pts = branch_feature_matrices_mp(args, sigma)
        K += mp.mpf("0.5") * branch_weight * (m.T * n + n.T * m)
        pieces.append((sigma, m, n, branch_weight, u_pts))
    K = sym(K)
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    _centers, R = trace_matrix(polys, local_args(args))
    F, trace_rank, trace_nullity, U = quotient_minimizer(K, R, args)

    plus_rows = None
    minus_rows = None
    row_us = []
    row_sigmas = []
    for sigma, m, n, weight, u_pts in pieces:
        scale = mp.sqrt(weight / 4)
        pblock = scale * (m + n) * F
        mblock = scale * (m - n) * F
        row_us.extend(u_pts)
        row_sigmas.extend([sigma for _ in u_pts])
        if plus_rows is None:
            plus_rows = pblock
            minus_rows = mblock
        else:
            plus_rows = stack_rows(plus_rows, pblock)
            minus_rows = stack_rows(minus_rows, mblock)
    P = sym(plus_rows.T * plus_rows)
    M = sym(minus_rows.T * minus_rows)
    C = plus_rows.T * minus_rows
    return {
        "traceRank": trace_rank,
        "traceNullity": trace_nullity,
        "P": P,
        "M": M,
        "C": C,
        "plusRows": plus_rows,
        "minusRows": minus_rows,
        "rowUs": row_us,
        "rowSigmas": row_sigmas,
        "R": R,
        "U": U,
    }


def stack_rows(a: mp.matrix, b: mp.matrix) -> mp.matrix:
    out = mp.matrix(a.rows + b.rows, a.cols)
    for i in range(a.rows):
        for j in range(a.cols):
            out[i, j] = a[i, j]
    for i in range(b.rows):
        for j in range(b.cols):
            out[a.rows + i, j] = b[i, j]
    return out


def hardy_targets(xi: XiTransform, omega: float, rmax: float | None):
    nodes = default_z_nodes()
    n = len(nodes)
    hp = [hardy_feature_scalar(xi, omega, z, "plus") for z in nodes]
    hm = [hardy_feature_scalar(xi, omega, z, "minus") for z in nodes]
    P = mp.matrix(n)
    M = mp.matrix(n)
    C = mp.matrix(n)
    for i, w in enumerate(nodes):
        for j, z in enumerate(nodes):
            P[i, j] = halfline_inner(hp[j], hp[i], z, w, rmax)
            M[i, j] = halfline_inner(hm[j], hm[i], z, w, rmax)
            C[i, j] = halfline_inner(hm[j], hp[i], z, w, rmax)
    return {"P": herm(P), "M": herm(M), "C": C, "joint": block2(herm(P), C, ctranspose(C), herm(M))}


def cholesky_star(a: mp.matrix) -> mp.matrix:
    L = mp.cholesky(herm(a))
    return ctranspose(L)


def random_orthonormal(rows: int, cols: int, rng: random.Random) -> mp.matrix:
    vecs = []
    for _ in range(cols):
        v = mp.matrix(rows, 1)
        for i in range(rows):
            v[i] = mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1))
        for q in vecs:
            v -= q * (ctranspose(q) * v)[0]
        norm = mp.sqrt((ctranspose(v) * v)[0])
        if norm == 0:
            return random_orthonormal(rows, cols, rng)
        vecs.append(v / norm)
    out = mp.matrix(rows, cols)
    for j, v in enumerate(vecs):
        for i in range(rows):
            out[i, j] = v[i]
    return out


def pinv_sqrt_real(P: mp.matrix, tol: mp.mpf):
    vals, vecs = mp.eigsy(sym(P), eigvals_only=False)
    maxv = max(abs(v) for v in vals) if len(vals) else mp.mpf("0")
    keep = [i for i, v in enumerate(vals) if v > tol * max(mp.mpf("1"), maxv)]
    V = columns(vecs, keep)
    D = mp.diag([1 / mp.sqrt(vals[i]) for i in keep])
    return V * D, len(keep), vals


def trace_coordinates_for_plus(model: dict, target: dict, Q: mp.matrix, args: argparse.Namespace):
    Pinv, rank, _vals = pinv_sqrt_real(model["P"], mp.mpf(args.rank_tol))
    B = cholesky_star(target["P"])
    return Pinv * Q * B, rank


def mapped_joint(model: dict, X: mp.matrix):
    P = herm(ctranspose(X) * model["P"] * X)
    M = herm(ctranspose(X) * model["M"] * X)
    C = ctranspose(X) * model["C"] * X
    return block2(P, C, ctranspose(C), M)


def block_residuals(mapped: mp.matrix, target: mp.matrix, n: int):
    return {
        "jointRelative": rel_frob(mapped, target),
        "plusRelative": rel_frob(mapped[:n, :n], target[:n, :n]),
        "crossRelative": rel_frob(mapped[:n, n:], target[:n, n:]),
        "minusRelative": rel_frob(mapped[n:, n:], target[n:, n:]),
        "jointFrobenius": frob_norm(mapped - target),
    }


def to_float(x) -> float:
    return float(mp.re(x))


def search(args: argparse.Namespace, model: dict, target: dict):
    rng = random.Random(args.seed)
    n = target["P"].rows
    _pinv, rank, _vals = pinv_sqrt_real(model["P"], mp.mpf(args.rank_tol))
    best = None
    best_x = None
    trials = max(1, args.random_trials)
    for trial in range(trials):
        if trial == 0:
            Q = mp.matrix(rank, n)
            for i in range(n):
                Q[i, i] = 1
        else:
            Q = random_orthonormal(rank, n, rng)
        X, _rank = trace_coordinates_for_plus(model, target, Q, args)
        mapped = mapped_joint(model, X)
        res = block_residuals(mapped, target["joint"], n)
        if best is None or res["jointRelative"] < best["jointRelative"]:
            best = {k: to_float(v) for k, v in res.items()}
            best["trial"] = trial
            best["traceCoordinateNormFrobenius"] = to_float(frob_norm(X))
            best_x = X
    return best, best_x, rank


def solve_regularized(a: mp.matrix, b: mp.matrix, reg: mp.mpf) -> mp.matrix:
    lhs = a.copy()
    for i in range(min(lhs.rows, lhs.cols)):
        lhs[i, i] += reg
    return mp.lu_solve(lhs, b)


def profile_target_vector(
    xi: XiTransform,
    omega: float,
    z: complex,
    rows_u: list[mp.mpf],
    branch: str,
) -> mp.matrix:
    amp = hardy_feature_scalar(xi, omega, z, branch)
    out = mp.matrix(len(rows_u), 1)
    for i, u in enumerate(rows_u):
        out[i] = amp * mp.e ** (1j * z * float(u))
    return out


def profile_fit_coordinates(args: argparse.Namespace, model: dict, xi: XiTransform, combined: bool):
    nodes = default_z_nodes()
    rank = model["traceRank"]
    X = mp.matrix(rank, len(nodes))
    P = sym(model["P"])
    M = sym(model["M"])
    Rplus = model["plusRows"]
    Rminus = model["minusRows"]
    reg = mp.mpf("1e-30")
    if combined:
        normal = P + M
    else:
        normal = P
    for j, z in enumerate(nodes):
        hp = profile_target_vector(xi, float(args.omega), z, model["rowUs"], "plus")
        rhs = Rplus.T * hp
        if combined:
            hm = profile_target_vector(xi, float(args.omega), z, model["rowUs"], "minus")
            rhs += Rminus.T * hm
        coeff = solve_regularized(normal, rhs, reg)
        for i in range(rank):
            X[i, j] = coeff[i]
    return X


def evaluate_x(model: dict, target: dict, X: mp.matrix):
    mapped = mapped_joint(model, X)
    return {k: to_float(v) for k, v in block_residuals(mapped, target["joint"], target["P"].rows).items()}


def preview_matrix(a: mp.matrix, max_rows: int = 5, max_cols: int = 5):
    rows = []
    for i in range(min(a.rows, max_rows)):
        row = []
        for j in range(min(a.cols, max_cols)):
            row.append({"real": float(mp.re(a[i, j])), "imag": float(mp.im(a[i, j]))})
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=5)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-24")
    parser.add_argument("--psd-tol", default="1e-24")
    parser.add_argument("--s-order", type=int, default=24)
    parser.add_argument("--u-order", type=int, default=50)
    parser.add_argument("--u-max", default="10")
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=1000)
    parser.add_argument("--hardy-rmax", type=float, default=40.0)
    parser.add_argument("--random-trials", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="klm_debranges_trace_map_constructor.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    model = finite_trace_model(args)
    full = hardy_targets(xi, float(args.omega), None)
    trunc = hardy_targets(xi, float(args.omega), args.hardy_rmax)

    best_full, x_full, plus_rank = search(args, model, full)
    best_trunc, x_trunc, _ = search(args, model, trunc)
    x_profile_plus = profile_fit_coordinates(args, model, xi, combined=False)
    x_profile_both = profile_fit_coordinates(args, model, xi, combined=True)
    profile_plus_full = evaluate_x(model, full, x_profile_plus)
    profile_both_full = evaluate_x(model, full, x_profile_both)
    profile_plus_trunc = evaluate_x(model, trunc, x_profile_plus)
    profile_both_trunc = evaluate_x(model, trunc, x_profile_both)
    threshold = 1e-6
    data = {
        "theoremName": "finite construction attempt for z -> x_z trace map",
        "omega": float(args.omega),
        "basis": args.basis,
        "constraints": args.constraints,
        "traceRank": model["traceRank"],
        "traceNullity": model["traceNullity"],
        "plusGramRank": plus_rank,
        "sOrder": args.s_order,
        "uOrder": args.u_order,
        "randomTrials": args.random_trials,
        "construction": (
            "X is chosen so X^*P_volterra X equals the Hardy plus Gram; "
            "remaining orthonormal freedom is searched for joint Gram matching."
        ),
        "fullHardyTarget": best_full,
        "truncatedHardyTarget": {**best_trunc, "hardyRmax": args.hardy_rmax},
        "profileFitCandidates": {
            "plusOnlyFull": profile_plus_full,
            "combinedPlusMinusFull": profile_both_full,
            "plusOnlyTruncated": profile_plus_trunc,
            "combinedPlusMinusTruncated": profile_both_trunc,
            "description": (
                "Profile fits solve least squares against sampled Hardy branch "
                "functions on the Volterra u-row grid, then check joint Grams."
            ),
        },
        "xFullPreview": preview_matrix(x_full),
        "xTruncatedPreview": preview_matrix(x_trunc),
        "exactTraceMapConstructed": best_full["jointRelative"] < threshold,
        "truncatedTraceMapConstructed": best_trunc["jointRelative"] < threshold,
        "diagnosis": (
            "Plus-Gram matching is not enough: the best finite trace coordinates "
            "still miss the Hardy minus/cross branch Grams.  The missing map "
            "must be derived from the actual Lambda_a trace action on the "
            "Hardy evaluation exponential."
            if best_full["jointRelative"] >= threshold
            else "Finite trace coordinates match the joint Hardy branch Grams."
        ),
        "nextProofTarget": (
            "Derive x_z(a)=Lambda_a applied to the primitive/evaluation "
            "exponential associated with h_z, then test those coordinates "
            "against the joint Hardy branch Grams."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Finite construction attempt for z -> x_z")
    print(f"  basis={args.basis} constraints={args.constraints} trace_rank={model['traceRank']} plus_rank={plus_rank}")
    print(f"  full joint rel residual: {best_full['jointRelative']:.6e}")
    print(
        f"    plus={best_full['plusRelative']:.3e} "
        f"cross={best_full['crossRelative']:.3e} minus={best_full['minusRelative']:.3e}"
    )
    print(f"  truncated R={args.hardy_rmax:g} joint rel residual: {best_trunc['jointRelative']:.6e}")
    print(
        f"    plus={best_trunc['plusRelative']:.3e} "
        f"cross={best_trunc['crossRelative']:.3e} minus={best_trunc['minusRelative']:.3e}"
    )
    print(f"  exact trace map constructed: {data['exactTraceMapConstructed']}")
    print(f"  truncated trace map constructed: {data['truncatedTraceMapConstructed']}")
    print(f"  profile plus-only full joint rel residual: {profile_plus_full['jointRelative']:.6e}")
    print(f"  profile combined full joint rel residual: {profile_both_full['jointRelative']:.6e}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
