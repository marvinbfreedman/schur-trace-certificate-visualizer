#!/usr/bin/env python3
import argparse
import math

from moment_transform_test import branch_matrix
from reduced_exact_finite import endpoint_ratios, pieces, psi

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("dilation_proof_diagnostics.py requires numpy; run with python") from exc


def quadrature(end, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return 0.5 * end * (nodes + 1.0), 0.5 * end * weights


def a_value(s, u, pcs):
    es = math.exp(s)
    eu = math.exp(u)
    total = 0.0
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * math.exp(beta * u - c * es * (eu - 1.0))
    return total


def log_psi(v, pcs):
    ev = math.exp(v)
    terms = []
    for coeff, beta, c in pcs:
        terms.append((1.0 if coeff > 0 else -1.0, math.log(abs(coeff)) + beta * v - c * ev))
    shift = max(logabs for _, logabs in terms)
    scaled = sum(sign * math.exp(logabs - shift) for sign, logabs in terms)
    if scaled <= 0.0:
        raise ValueError(f"psi is not positive at v={v:g}")
    return shift + math.log(scaled)


def center_layer_matrix(kind, center, interval_end, order):
    pcs = pieces(kind)
    pts, wts = quadrature(interval_end, order)
    root = np.sqrt(wts)
    mat = np.zeros((order, order), dtype=float)
    for i, s in enumerate(pts):
        log_ps = log_psi(float(s), pcs)
        for j, t in enumerate(pts[: i + 1]):
            value = 0.0
            if s + t <= 2.0 * center:
                left = center + 0.5 * (s - t)
                right = center - 0.5 * (s - t)
                if left >= 0.0 and right >= 0.0:
                    value = math.exp(
                        log_psi(float(left), pcs)
                        + log_psi(float(right), pcs)
                        - log_ps
                        - log_psi(float(t), pcs)
                    )
            mat[i, j] = mat[j, i] = root[i] * value * root[j]
    return mat


def branch_second_derivative_matrix(kind, exponent, interval_end, s_order, u_end, u_order):
    pcs = pieces(kind)
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root = np.sqrt(s_wts)
    mat = np.zeros((s_order, s_order), dtype=float)
    for u, w in zip(u_pts, u_wts):
        values = np.array(
            [
                a_value(float(s), float(u), pcs) * math.exp(exponent * float(s + u))
                for s in s_pts
            ],
            dtype=float,
        )
        x = s_pts + u
        layer = (
            np.outer(x * x * values, values)
            + 2.0 * np.outer(x * values, x * values)
            + np.outer(values, x * x * values)
        )
        mat += float(w) * root[:, None] * layer * root[None, :]
    return mat


def paired_second_derivative_matrix(kind, omega, lam, interval_end, s_order, u_end, u_order):
    if omega == 0.0:
        return branch_second_derivative_matrix(
            kind, lam, interval_end, s_order, u_end, u_order
        )
    half = 0.5 * omega
    return 0.5 * (
        branch_second_derivative_matrix(
            kind, lam + half, interval_end, s_order, u_end, u_order
        )
        + branch_second_derivative_matrix(
            kind, lam - half, interval_end, s_order, u_end, u_order
        )
    )


def branch_first_derivative_matrix(kind, exponent, interval_end, s_order, u_end, u_order):
    if abs(exponent) < 1e-15:
        omega = 0.0
        sigma = 0
    else:
        omega = 2.0 * abs(exponent)
        sigma = 1 if exponent > 0 else -1
    m, n = branch_matrix(kind, omega, sigma, interval_end, s_order, u_end, u_order)
    return 0.5 * (m.T @ n + n.T @ m)


def paired_first_derivative_matrix(kind, omega, lam, interval_end, s_order, u_end, u_order):
    if omega == 0.0:
        return branch_first_derivative_matrix(
            kind, lam, interval_end, s_order, u_end, u_order
        )
    half = 0.5 * omega
    return 0.5 * (
        branch_first_derivative_matrix(kind, lam + half, interval_end, s_order, u_end, u_order)
        + branch_first_derivative_matrix(kind, lam - half, interval_end, s_order, u_end, u_order)
    )


def eig_report(name, mat, tol):
    vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
    print(
        f"  {name}: min={vals[0]:.12e} max={vals[-1]:.12e} "
        f"neg={(vals < -tol).sum()}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=60)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=200)
    parser.add_argument("--centers", default="0.1,0.2,0.5,1,2,4")
    parser.add_argument("--lambdas", default="0,0.05,0.1,0.2,0.5")
    parser.add_argument("--tol", type=float, default=1e-10)
    args = parser.parse_args()

    print(
        f"dilation proof diagnostics kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    print("center-layer test for x = u + (s+t)/2")
    for text in args.centers.split(","):
        center = float(text)
        eig_report(
            f"center={center:g}",
            center_layer_matrix(args.kind, center, args.L, args.s_order),
            args.tol,
        )

    print("paired first derivative D_lambda = N'_lambda / 2")
    for text in args.lambdas.split(","):
        lam = float(text)
        eig_report(
            f"lambda={lam:g}",
            paired_first_derivative_matrix(
                args.kind, args.omega, lam, args.L, args.s_order, args.umax, args.u_order
            ),
            args.tol,
        )

    print("paired second derivative N''_lambda")
    for text in args.lambdas.split(","):
        lam = float(text)
        eig_report(
            f"lambda={lam:g}",
            paired_second_derivative_matrix(
                args.kind, args.omega, lam, args.L, args.s_order, args.umax, args.u_order
            ),
            args.tol,
        )


if __name__ == "__main__":
    main()
