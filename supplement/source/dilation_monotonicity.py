#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import endpoint_ratios, pieces

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("dilation_monotonicity.py requires numpy; run with python, not python3") from exc


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


def branch_transform_matrix(pcs, omega, lam, sigma, interval_end, s_order, u_end, u_order):
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root_s = np.sqrt(s_wts)
    bmat = np.zeros((u_order, s_order), dtype=float)
    for i, u in enumerate(u_pts):
        root_u = math.sqrt(float(u_wts[i]))
        for j, s in enumerate(s_pts):
            branch = math.exp((0.5 * sigma * omega + lam) * (s + u))
            bmat[i, j] = root_u * a_value(float(s), float(u), pcs) * branch * root_s[j]
    return bmat, s_pts, u_pts


def norm_matrix(kind, omega, lam, interval_end, s_order, u_end, u_order):
    pcs = pieces(kind)
    mat = np.zeros((s_order, s_order), dtype=float)
    sigmas = [0] if omega == 0.0 else [-1, 1]
    branch_weights = [1.0] if omega == 0.0 else [0.5, 0.5]
    for sigma, bw in zip(sigmas, branch_weights):
        bmat, _, _ = branch_transform_matrix(
            pcs, omega, lam, sigma, interval_end, s_order, u_end, u_order
        )
        mat += bw * (bmat.T @ bmat)
    return mat


def moment_matrix(kind, omega, interval_end, s_order, u_end, u_order):
    pcs = pieces(kind)
    mat = np.zeros((s_order, s_order), dtype=float)
    sigmas = [0] if omega == 0.0 else [-1, 1]
    branch_weights = [1.0] if omega == 0.0 else [0.5, 0.5]
    for sigma, bw in zip(sigmas, branch_weights):
        mmat, s_pts, u_pts = branch_transform_matrix(
            pcs, omega, 0.0, sigma, interval_end, s_order, u_end, u_order
        )
        factor = np.add.outer(u_pts, s_pts)
        nmat = factor * mmat
        mat += 0.5 * bw * (mmat.T @ nmat + nmat.T @ mmat)
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=70)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=220)
    parser.add_argument("--lambdas", default="0,0.001,0.005,0.01,0.02,0.05")
    parser.add_argument("--derivative-step", type=float, default=1e-5)
    args = parser.parse_args()

    lambdas = [float(text) for text in args.lambdas.split(",")]
    mats = [
        norm_matrix(args.kind, args.omega, lam, args.L, args.s_order, args.umax, args.u_order)
        for lam in lambdas
    ]
    print(
        f"dilation monotonicity kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    for lam, mat in zip(lambdas, mats):
        vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
        print(f"  lambda={lam:g}: norm-kernel min={vals[0]:.12e} max={vals[-1]:.12e}")
    for (lam0, mat0), (lam1, mat1) in zip(zip(lambdas, mats), zip(lambdas[1:], mats[1:])):
        diff = mat1 - mat0
        vals = np.linalg.eigvalsh((diff + diff.T) / 2.0)
        print(
            f"  diff {lam1:g}-{lam0:g}: min={vals[0]:.12e} "
            f"max={vals[-1]:.12e} neg={(vals < -1e-12).sum()}"
        )
    qmat = moment_matrix(args.kind, args.omega, args.L, args.s_order, args.umax, args.u_order)
    qvals = np.linalg.eigvalsh((qmat + qmat.T) / 2.0)
    h = args.derivative_step
    fd = (
        norm_matrix(args.kind, args.omega, h, args.L, args.s_order, args.umax, args.u_order)
        - mats[0]
    ) / h
    defect = fd - 2.0 * qmat
    scale = max(1.0, float(np.max(np.abs(2.0 * qmat))))
    print(f"  moment derivative Q_min={qvals[0]:.12e} Q_max={qvals[-1]:.12e}")
    print(
        f"  derivative check: ||(N_h-N_0)/h - 2Q||_max="
        f"{np.max(np.abs(defect)):.12e} rel={np.max(np.abs(defect)) / scale:.12e}"
    )


if __name__ == "__main__":
    main()
