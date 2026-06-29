#!/usr/bin/env python3
import argparse
import math

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("positive_branch_perturbation.py requires numpy; run with python") from exc


PI = math.pi


def quadrature(end, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return 0.5 * end * (nodes + 1.0), 0.5 * end * weights


def psi_mode(n, v):
    c = PI * n * n
    ev = math.exp(v)
    return (4.0 * c * c * math.exp(2.25 * v) - 6.0 * c * math.exp(1.25 * v)) * math.exp(
        -c * ev
    )


def branch_matrix(mode_weights, exponent, interval_end, s_order, u_end, u_order):
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root = np.sqrt(s_wts)
    mat = np.zeros((s_order, s_order), dtype=float)
    for u, w in zip(u_pts, u_wts):
        values = []
        for s in s_pts:
            v = float(s + u)
            total = sum(weight * psi_mode(mode, v) for mode, weight in mode_weights.items())
            values.append(math.exp(exponent * v) * total)
        values = np.array(values, dtype=float)
        x = s_pts + u
        layer = 0.5 * (np.outer(x * values, values) + np.outer(values, x * values))
        mat += float(w) * root[:, None] * layer * root[None, :]
    return mat


def min_eig(mode_weights, exponent, interval_end, s_order, u_end, u_order):
    mat = branch_matrix(mode_weights, exponent, interval_end, s_order, u_end, u_order)
    vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
    return vals[0], vals[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=float, default=0.245)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=50)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=180)
    parser.add_argument("--rmax", type=float, default=20.0)
    parser.add_argument("--tol", type=float, default=1e-10)
    args = parser.parse_args()

    print(
        f"positive branch perturbation a={args.a:g} L={args.L:g} "
        f"s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    for weights in ({1: 1.0}, {2: 1.0}, {1: 1.0, 2: 1.0}, {1: 1.0, 2: 10.0}):
        low, high = min_eig(weights, args.a, args.L, args.s_order, args.umax, args.u_order)
        print(f"  weights={weights}: min={low:.12e} max={high:.12e}")

    lo = 0.0
    hi = args.rmax
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        low, _ = min_eig({1: 1.0, 2: mid}, args.a, args.L, args.s_order, args.umax, args.u_order)
        if low >= -args.tol:
            lo = mid
        else:
            hi = mid
    lo_min, _ = min_eig({1: 1.0, 2: lo}, args.a, args.L, args.s_order, args.umax, args.u_order)
    hi_min, _ = min_eig({1: 1.0, 2: hi}, args.a, args.L, args.s_order, args.umax, args.u_order)
    print(
        f"  visible-negativity threshold for psi_1 + r psi_2: "
        f"r in [{lo:.8g}, {hi:.8g}]"
    )
    print(f"    min at lower={lo_min:.12e}, upper={hi_min:.12e}")


if __name__ == "__main__":
    main()
