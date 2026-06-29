#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import endpoint_ratios, pieces

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("moment_transform_test.py requires numpy; run with python, not python3") from exc


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


def branch_matrix(kind, omega, sigma, interval_end, s_order, u_end, u_order):
    pcs = pieces(kind)
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root_s = np.sqrt(s_wts)
    root_u = np.sqrt(u_wts)
    m = np.zeros((u_order, s_order), dtype=float)
    n = np.zeros_like(m)
    for i, u in enumerate(u_pts):
        for j, s in enumerate(s_pts):
            branch = math.exp(0.5 * sigma * omega * (s + u))
            value = root_u[i] * a_value(float(s), float(u), pcs) * branch * root_s[j]
            m[i, j] = value
            n[i, j] = float(s + u) * value
    return m, n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=80)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=240)
    args = parser.parse_args()

    total_q = np.zeros((args.s_order, args.s_order), dtype=float)
    total_h = np.zeros_like(total_q)
    print(
        f"moment transform kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    sigmas = [0] if args.omega == 0.0 else [-1, 1]
    weights = [1.0] if args.omega == 0.0 else [0.5, 0.5]
    for sigma, branch_weight in zip(sigmas, weights):
        m, n = branch_matrix(
            args.kind, args.omega, sigma, args.L, args.s_order, args.umax, args.u_order
        )
        q = 0.5 * branch_weight * (m.T @ n + n.T @ m)
        h = branch_weight * (m.T @ m)
        total_q += q
        total_h += h
        q_vals = np.linalg.eigvalsh((q + q.T) / 2.0)
        h_vals = np.linalg.eigvalsh((h + h.T) / 2.0)
        print(
            f"  sigma={sigma}: Q_min={q_vals[0]:.12e} Q_max={q_vals[-1]:.12e} "
            f"H_min={h_vals[0]:.12e} H_rank={(h_vals > 1e-12).sum()}"
        )
    q_vals = np.linalg.eigvalsh((total_q + total_q.T) / 2.0)
    h_vals = np.linalg.eigvalsh((total_h + total_h.T) / 2.0)
    print(
        f"  total: Q_min={q_vals[0]:.12e} Q_max={q_vals[-1]:.12e} "
        f"H_min={h_vals[0]:.12e} H_rank={(h_vals > 1e-12).sum()}"
    )
    print("  Q low=", " ".join(f"{v:.4e}" for v in q_vals[:10]))


if __name__ == "__main__":
    main()
