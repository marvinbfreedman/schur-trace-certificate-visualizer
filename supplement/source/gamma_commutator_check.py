#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("gamma_commutator_check.py requires numpy; run with python") from exc


def base_kernel(p, c, v):
    return math.exp(p * v - c * math.exp(v))


def base_kernel_derivative(p, c, v):
    return (p - c * math.exp(v)) * base_kernel(p, c, v)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=1.495)
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--L", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=240)
    args = parser.parse_args()

    pts, wts = quadrature(args.L, args.order)
    g = np.exp(-0.7 * pts) * (1.0 + 0.4 * pts)
    dg = -0.7 * np.exp(-0.7 * pts) * (1.0 + 0.4 * pts) + 0.4 * np.exp(-0.7 * pts)
    print(
        f"gamma commutator check p={args.p:g} c={args.c:g} "
        f"L={args.L:g} order={args.order}"
    )
    for r in (0.0, 0.1, 1.0, 3.0):
        dh = sum(
            w * base_kernel_derivative(args.p, args.c, float(r + s)) * gv
            for s, w, gv in zip(pts, wts, g)
        )
        hd = sum(
            w * base_kernel(args.p, args.c, float(r + s)) * dgv
            for s, w, dgv in zip(pts, wts, dg)
        )
        boundary = -base_kernel(args.p, args.c, r) * 1.0
        print(
            f"  r={r:g}: (D H + H D)g={dh + hd:.12e} "
            f"boundary={boundary:.12e} defect={dh + hd - boundary:.12e}"
        )


if __name__ == "__main__":
    main()
