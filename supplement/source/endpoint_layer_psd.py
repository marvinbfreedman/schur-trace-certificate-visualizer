#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_layer_psd.py requires numpy; run with python") from exc


def k_value(p, c, x):
    return 2.0 * c * math.exp(p * x) * (2.0 * c * math.exp(x) - 3.0) * math.exp(
        -c * math.exp(x)
    )


def layer_matrix(p, c, interval_end, s_order, alpha, u_end, u_order):
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root = np.sqrt(s_wts)
    threshold = 0.5 * (s_pts[:, None] + s_pts[None, :])
    mat = np.zeros((s_order, s_order), dtype=float)
    for u, w in zip(u_pts, u_wts):
        values = np.array([k_value(p, c, float(u + s)) for s in s_pts])
        active = (float(u) + threshold) >= alpha
        mat += float(w) * active * np.outer(values, values)
    return root[:, None] * mat * root[None, :]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=1.25)
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--L", type=float, default=12.0)
    parser.add_argument("--s-order", type=int, default=50)
    parser.add_argument("--umax", type=float, default=18.0)
    parser.add_argument("--u-order", type=int, default=180)
    parser.add_argument("--alphas", type=float, nargs="*", default=None)
    parser.add_argument("--alpha-max", type=float, default=10.0)
    parser.add_argument("--alpha-count", type=int, default=21)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    if args.alphas is None:
        alphas = [
            args.alpha_max * i / (args.alpha_count - 1)
            for i in range(args.alpha_count)
        ]
    else:
        alphas = args.alphas
    print(
        f"endpoint layer PSD p={args.p:g} L={args.L:g} s_order={args.s_order} "
        f"u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    worst = (float("inf"), None)
    for alpha in alphas:
        mat = layer_matrix(
            args.p, args.c, args.L, args.s_order, alpha, args.umax, args.u_order
        )
        vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
        if vals[0] < worst[0]:
            worst = (float(vals[0]), alpha)
        print(
            f"  alpha={alpha:.6g} min={vals[0]: .12e} "
            f"max={vals[-1]: .12e} neg={(vals < -args.tol).sum()}"
        )
    print(f"worst min={worst[0]:.12e} at alpha={worst[1]:.6g}")


if __name__ == "__main__":
    main()
