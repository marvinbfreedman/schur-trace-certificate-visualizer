#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_log_transport.py requires numpy; run with python") from exc


def h_value(c, u):
    eu = math.exp(u)
    return math.exp(1.25 * u) * (c * eu - 1.5) * math.exp(-c * eu)


def p_kernel(c, t, s, r_end, r_order):
    pts, wts = quadrature(r_end, r_order)
    total = 0.0
    for r, w in zip(pts, wts):
        total += (
            float(w)
            * (float(r) + 0.5 * (t + s))
            * h_value(c, t + float(r))
            * h_value(c, s + float(r))
        )
    return total


def transport_source(c, t, s):
    return -0.5 * (t + s) * h_value(c, t) * h_value(c, s)


def finite_transport_derivative(c, t, s, eps, r_end, r_order):
    return (
        p_kernel(c, t + eps, s + eps, r_end, r_order)
        - p_kernel(c, t, s, r_end, r_order)
    ) / eps


def eigvals(mat):
    return np.linalg.eigvalsh((mat + mat.T) / 2.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--r-order", type=int, default=260)
    parser.add_argument("--eps", type=float, default=1e-5)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    pts, wts = quadrature(args.T, args.order)
    root = np.sqrt(wts)
    pmat = np.zeros((args.order, args.order), dtype=float)
    smat = np.zeros_like(pmat)
    for i, t in enumerate(pts):
        for j, s in enumerate(pts[: i + 1]):
            p = p_kernel(args.c, float(t), float(s), args.rmax, args.r_order)
            src = -transport_source(args.c, float(t), float(s))
            pmat[i, j] = pmat[j, i] = root[i] * p * root[j]
            smat[i, j] = smat[j, i] = root[i] * src * root[j]

    pvals = eigvals(pmat)
    svals = eigvals(smat)
    print(
        f"endpoint log transport T=[0,{args.T:g}] order={args.order} "
        f"r=[0,{args.rmax:g}] r_order={args.r_order}"
    )
    print(
        f"  P kernel min={pvals[0]:.12e} max={pvals[-1]:.12e} "
        f"neg={(pvals < -args.tol).sum()}"
    )
    print(
        f"  -transport source min={svals[0]:.12e} max={svals[-1]:.12e} "
        f"neg={(svals < -args.tol).sum()}"
    )
    print("  transport checks:")
    for t, s in ((0.0, 0.0), (0.0, 0.7), (0.4, 1.1), (1.0, 1.8)):
        fd = finite_transport_derivative(
            args.c, t, s, args.eps, args.rmax, args.r_order
        )
        exact = transport_source(args.c, t, s)
        print(
            f"    t={t:g} s={s:g}: finite={fd:.12e} "
            f"exact={exact:.12e} defect={fd - exact:.12e}"
        )


if __name__ == "__main__":
    main()
