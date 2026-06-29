#!/usr/bin/env python3
"""Construct a finite-dimensional quotient factorization certificate.

For a Galerkin matrix K for the reduced Volterra form and a sampled endpoint
trace matrix

  (R f)(a_j) = Lambda_{a_j}(f),

this script checks the exact finite-dimensional condition behind the proposed
quotient theorem:

  K = P - R^* D R,        P >= 0, D >= 0.

Equivalently, K must be nonnegative on ker R and the coupling from ker R to
the row-space of R must satisfy the usual range condition.  When those hold,
we construct D explicitly.  In the continuum proof, this is the model for

  Q_Psi(f) = ||Gf||^2 - ||S Rf||^2.

The implementation is pure mpmath because numpy is not always available in
this workspace.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at  # noqa: E402
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


def poly_add(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf("0") for _ in range(n)]
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def poly_scale(a, c):
    return trim([c * x for x in a])


def poly_mul(a, b):
    out = [mp.mpf("0") for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def trim(a):
    while len(a) > 1 and abs(a[-1]) == 0:
        a.pop()
    return a


def shifted_legendre_polys(order, length):
    z = [mp.mpf("-1"), mp.mpf("2") / length]
    polys = [[mp.mpf("1")]]
    if order == 1:
        return [poly_scale(polys[0], 1 / mp.sqrt(length))]
    polys.append(z)
    for n in range(1, order - 1):
        top = poly_add(
            poly_scale(poly_mul(z, polys[n]), 2 * n + 1),
            poly_scale(polys[n - 1], -n),
        )
        polys.append(poly_scale(top, mp.mpf(1) / (n + 1)))
    out = []
    for n, poly in enumerate(polys):
        out.append(poly_scale(poly, mp.sqrt(mp.mpf(2 * n + 1) / length)))
    return out


def poly_value(poly, x):
    total = mp.mpf("0")
    for coeff in reversed(poly):
        total = total * x + coeff
    return total


def poly_derivative_value(poly, x, deriv):
    total = mp.mpf("0")
    for power in range(deriv, len(poly)):
        coeff = poly[power] * mp.factorial(power) / mp.factorial(power - deriv)
        total += coeff * x ** (power - deriv)
    return total


def legendre_quadrature(length, order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    pts = []
    wts = []
    for node, weight in zip(nodes, weights):
        pts.append(mp.mpf("0.5") * length * (node + 1))
        wts.append(mp.mpf("0.5") * length * weight)
    return pts, wts


def psi_value(v, pieces):
    ev = mp.e**v
    return mp.fsum(coeff * mp.e ** (beta * v - c * ev) for coeff, beta, c in pieces)


def endpoint_ratios(v, pieces):
    denom = psi_value(v, pieces)
    ev = mp.e**v
    return [
        (coeff * mp.e ** (beta * v - c * ev) / denom, beta, c)
        for coeff, beta, c in pieces
    ]


def laguerre_integral(alpha, p, center, omega, nodes, weights):
    total = mp.mpf("0")
    for node, weight in zip(nodes, weights):
        q = 1 + node / alpha
        u = center + mp.log(q)
        total += weight * u * mp.cosh(omega * u) * q**p
    return total / alpha


def kernel_red(s, t, omega, pieces, lag_nodes, lag_weights):
    center = (s + t) / 2
    left = endpoint_ratios(s, pieces)
    right = endpoint_ratios(t, pieces)
    es = mp.e**s
    et = mp.e**t
    total = mp.mpf("0")
    for ratio_i, beta_i, c_i in left:
        for ratio_j, beta_j, c_j in right:
            alpha = c_i * es + c_j * et
            p = beta_i + beta_j - 1
            total += ratio_i * ratio_j * laguerre_integral(
                alpha, p, center, omega, lag_nodes, lag_weights
            )
    return total


def endpoint_b_vw(s, r, c):
    x = mp.e**r
    tau = x - 1
    lam = c * mp.e**s
    denom = lam - c
    exp_lam = mp.e ** (-lam * tau)
    exp_c = mp.e ** (-c * tau)
    v = ((1 - lam * x) * exp_lam - (1 - c * x) * exp_c) / denom
    w = s * (1 - lam * x) * exp_lam / denom
    return v, w


def segments(rmax):
    pts = [mp.mpf("0"), mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")]
    z = mp.mpf("0.4")
    while z < rmax:
        pts.append(z)
        z *= 2
    pts.append(rmax)
    return pts


def endpoint_b_quadrature(rmax, order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    out_nodes = []
    out_weights = []
    for a, b in zip(segments(rmax)[:-1], segments(rmax)[1:]):
        mid = (a + b) / 2
        half = (b - a) / 2
        for node, weight in zip(nodes, weights):
            out_nodes.append(mid + half * node)
            out_weights.append(half * weight)
    return out_nodes, out_weights


def kernel_endpoint_b(s, t, c, r_nodes, r_weights):
    total = mp.mpf("0")
    for r, weight in zip(r_nodes, r_weights):
        vs, ws = endpoint_b_vw(s, r, c)
        vt, wt = endpoint_b_vw(t, r, c)
        total += weight * mp.e ** (mp.mpf("2.5") * r) * (
            r * vs * vt + mp.mpf("0.5") * (vs * wt + ws * vt)
        )
    return total


def gram_matrix(args, omega, length):
    kind = args.kind
    basis_order = args.basis
    quad_order = args.quad
    pieces = pieces_for(kind)
    polys = shifted_legendre_polys(basis_order, length)
    pts, wts = legendre_quadrature(length, quad_order)
    lag_nodes, lag_weights = mp.gauss_quadrature(args.laguerre, "laguerre")
    r_nodes, r_weights = endpoint_b_quadrature(
        mp.mpf(args.endpoint_kernel_rmax), args.endpoint_kernel_order
    )

    values = mp.matrix(quad_order, basis_order)
    for i, x in enumerate(pts):
        for j, poly in enumerate(polys):
            values[i, j] = poly_value(poly, x)

    out = mp.matrix(basis_order)
    for a, s in enumerate(pts):
        for b, t in enumerate(pts):
            if args.model == "full":
                kval = kernel_red(s, t, omega, pieces, lag_nodes, lag_weights)
            elif args.model == "endpoint_b":
                kval = kernel_endpoint_b(s, t, mp.pi, r_nodes, r_weights)
            else:
                raise ValueError(args.model)
            weight = wts[a] * wts[b] * kval
            for i in range(basis_order):
                vi = values[a, i]
                for j in range(basis_order):
                    out[i, j] += weight * vi * values[b, j]
    return (out + out.T) / 2, polys


def sampled_centers(lo, hi, count):
    if count == 1:
        return [(lo + hi) / 2]
    step = (hi - lo) / (count - 1)
    return [lo + i * step for i in range(count)]


def signed_defect(center, jet_order, rmax, endpoint_order, tol):
    vals, neg, vec = defect_at(center, jet_order, rmax, endpoint_order, tol)
    if vec[0] > 0:
        vec = [-x for x in vec]
    return vals, neg, vec


def trace_matrix(polys, args):
    centers = sampled_centers(
        mp.mpf(args.constraint_min), mp.mpf(args.constraint_max), args.constraints
    )
    rows = []
    for center in centers:
        vals, neg, vec = signed_defect(
            center,
            args.jet_order,
            mp.mpf(args.endpoint_rmax),
            args.endpoint_order,
            mp.mpf(args.endpoint_tol),
        )
        row = []
        for poly in polys:
            value = mp.mpf("0")
            for k, coeff in enumerate(vec):
                value += coeff * poly_derivative_value(poly, center, k) / mp.factorial(k)
            row.append(value)
        rows.append(row)
    return centers, mp.matrix(rows)


def columns(mat, indices):
    out = mp.matrix(mat.rows, len(indices))
    for j, idx in enumerate(indices):
        for i in range(mat.rows):
            out[i, j] = mat[i, idx]
    return out


def diag(values):
    out = mp.matrix(len(values))
    for i, value in enumerate(values):
        out[i, i] = value
    return out


def frob_norm(mat):
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def op_norm(mat):
    if mat.rows == 0 or mat.cols == 0:
        return mp.mpf("0")
    vals = mp.eigsy((mat.T * mat + (mat.T * mat).T) / 2, eigvals_only=True)
    return mp.sqrt(max(mp.mpf("0"), vals[-1]))


def minmax_eigs(mat):
    vals = mp.eigsy((mat + mat.T) / 2, eigvals_only=True)
    return vals, vals[0], vals[-1]


def positive_part_inverse(mat, tol):
    vals, vecs = mp.eigsy((mat + mat.T) / 2, eigvals_only=False)
    keep = [i for i, val in enumerate(vals) if val > tol]
    zero = [i for i, val in enumerate(vals) if val <= tol]
    inv = mp.matrix(mat.rows)
    if keep:
        vkeep = columns(vecs, keep)
        inv = vkeep * diag([1 / vals[i] for i in keep]) * vkeep.T
    return vals, vecs, keep, zero, inv


def max_eig_or_zero(mat):
    if mat.rows == 0:
        return mp.mpf("0")
    vals = mp.eigsy((mat + mat.T) / 2, eigvals_only=True)
    return max(mp.mpf("0"), vals[-1])


def quotient_certificate(K, R, args):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)

    A = N.T * K * N
    B = N.T * K * U
    C = U.T * K * U

    avals, avecs, a_keep, a_zero, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    range_resid = mp.mpf("0")
    if a_zero and B.cols:
        zvecs = columns(avecs, a_zero)
        range_resid = frob_norm(zvecs.T * B)

    gamma2 = B.T * Aplus * B if B.cols else mp.matrix(0)
    H = gamma2 - C if B.cols else -C
    alpha = max_eig_or_zero(H) + mp.mpf(args.margin)
    M = alpha * mp.eye(U.cols)

    P_block = mp.matrix(K.rows)
    # In the [N,U] basis, P has lower block C + M and K has lower block C.
    basis = N
    if U.cols:
        basis = N.copy()
        basis = mp.matrix([[basis[i, j] for j in range(basis.cols)] + [U[i, j] for j in range(U.cols)] for i in range(K.rows)])
    P_in_basis = mp.matrix(K.rows)
    nn = N.cols
    uu = U.cols
    for i in range(nn):
        for j in range(nn):
            P_in_basis[i, j] = A[i, j]
    for i in range(nn):
        for j in range(uu):
            P_in_basis[i, nn + j] = B[i, j]
            P_in_basis[nn + j, i] = B[i, j]
    for i in range(uu):
        for j in range(uu):
            P_in_basis[nn + i, nn + j] = C[i, j] + M[i, j]
    P = basis * P_in_basis * basis.T

    # Construct D on sampled trace space.  If R U = W Sigma, then choosing
    # D = W Sigma^{-1} M Sigma^{-1} W^T gives U^T R^T D R U = M.
    D = mp.matrix(R.rows)
    if U.cols:
        sigmas = [mp.sqrt(rvals[i]) for i in u_idx]
        RU = R * U
        W = mp.matrix(R.rows, U.cols)
        for j, sigma in enumerate(sigmas):
            for i in range(R.rows):
                W[i, j] = RU[i, j] / sigma
        D = W * diag([alpha / (sigma * sigma) for sigma in sigmas]) * W.T

    P_from_trace = K + R.T * D * R
    pvals, pmin, pmax = minmax_eigs(P)
    p2vals, p2min, p2max = minmax_eigs(P_from_trace)
    dvals, dmin, dmax = minmax_eigs(D) if D.rows else ([], mp.mpf("0"), mp.mpf("0"))
    gvals, gmin, gmax = minmax_eigs(gamma2) if gamma2.rows else ([], mp.mpf("0"), mp.mpf("0"))
    hvals, hmin, hmax = minmax_eigs(H) if H.rows else ([], mp.mpf("0"), mp.mpf("0"))
    kvals, kmin, kmax = minmax_eigs(K)
    avals_full, amin, amax = minmax_eigs(A) if A.rows else ([], mp.mpf("0"), mp.mpf("0"))
    schur = C + M - B.T * Aplus * B if U.cols else mp.matrix(0)
    svals, smin, smax = minmax_eigs(schur) if U.cols else ([], mp.mpf("0"), mp.mpf("0"))
    trace_error = frob_norm(P - P_from_trace)
    b_norm = op_norm(B)
    normalized_range_resid = range_resid / b_norm if b_norm else mp.mpf("0")

    return {
        "rvals": rvals,
        "rank_tol": rank_tol,
        "nullity": N.cols,
        "rank": U.cols,
        "kmin": kmin,
        "kmax": kmax,
        "amin": amin,
        "amax": amax,
        "b_norm": b_norm,
        "range_resid": range_resid,
        "normalized_range_resid": normalized_range_resid,
        "gamma2_min": gmin,
        "gamma2_max": gmax,
        "h_min": hmin,
        "h_max": hmax,
        "alpha": alpha,
        "schur_min": smin,
        "schur_max": smax,
        "pmin": pmin,
        "pmax": pmax,
        "p2min": p2min,
        "p2max": p2max,
        "dmin": dmin,
        "dmax": dmax,
        "trace_error": trace_error,
        "low_k": kvals[: min(6, len(kvals))],
        "low_A": avals_full[: min(6, len(avals_full))],
        "low_gamma2": gvals[: min(6, len(gvals))],
        "low_H": hvals[: min(6, len(hvals))],
        "low_P": p2vals[: min(6, len(p2vals))],
        "low_D": dvals[: min(6, len(dvals))] if D.rows else [],
    }


def fmt(value, digits=14):
    return mp.nstr(value, digits)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="full")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=12)
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=36)
    parser.add_argument("--endpoint-kernel-order", type=int, default=24)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--constraints", type=int, default=5)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=35)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-30")
    parser.add_argument("--psd-tol", default="1e-30")
    parser.add_argument("--margin", default="1e-28")
    parser.add_argument("--dps", type=int, default=60)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    length = mp.mpf(args.L)
    omega = mp.mpf(args.omega)

    print(
        f"quotient factorization model={args.model} kind={args.kind} "
        f"omega={omega} L={length} "
        f"basis={args.basis} quad={args.quad} laguerre={args.laguerre}"
    )
    print(
        f"  constraints={args.constraints} interval=[{args.constraint_min},"
        f"{args.constraint_max}] dps={args.dps}"
    )
    K, polys = gram_matrix(args, omega, length)
    centers, R = trace_matrix(polys, args)
    cert = quotient_certificate(K, R, args)

    print("  centers=" + ", ".join(fmt(c, 8) for c in centers))
    print(
        f"  trace rank={cert['rank']} nullity={cert['nullity']} "
        f"rank_tol={fmt(cert['rank_tol'])}"
    )
    print(f"  K eig min={fmt(cert['kmin'])} max={fmt(cert['kmax'])}")
    print(f"  K low=" + " ".join(fmt(v, 8) for v in cert["low_k"]))
    print(
        f"  K|kerR eig min={fmt(cert['amin'])} max={fmt(cert['amax'])} "
        f"range_resid={fmt(cert['range_resid'])}"
    )
    print(
        f"  cross ||B||={fmt(cert['b_norm'])} "
        f"normalized_range_resid={fmt(cert['normalized_range_resid'])}"
    )
    print(f"  K|kerR low=" + " ".join(fmt(v, 8) for v in cert["low_A"]))
    print(
        f"  Douglas Gamma^*Gamma eig min={fmt(cert['gamma2_min'])} "
        f"max={fmt(cert['gamma2_max'])}"
    )
    print(
        "  Douglas low="
        + " ".join(fmt(v, 8) for v in cert["low_gamma2"])
    )
    print(
        f"  Schur defect Gamma^*Gamma-C eig min={fmt(cert['h_min'])} "
        f"max={fmt(cert['h_max'])}"
    )
    print(
        "  Schur defect low="
        + " ".join(fmt(v, 8) for v in cert["low_H"])
    )
    print(
        f"  constructed alpha={fmt(cert['alpha'])} "
        f"Schur_min={fmt(cert['schur_min'])}"
    )
    print(
        f"  D eig min={fmt(cert['dmin'])} max={fmt(cert['dmax'])} "
        f"low=" + " ".join(fmt(v, 8) for v in cert["low_D"])
    )
    print(
        f"  P=K+R^*DR eig min={fmt(cert['p2min'])} max={fmt(cert['p2max'])} "
        f"trace_error={fmt(cert['trace_error'])}"
    )
    print(f"  P low=" + " ".join(fmt(v, 8) for v in cert["low_P"]))


if __name__ == "__main__":
    main()
