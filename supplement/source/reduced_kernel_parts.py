#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import endpoint_ratios, pieces

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("reduced_kernel_parts.py requires numpy; run with python, not python3") from exc


def a_value(s, u, pcs, omega=0.0, sigma=0):
    total = 0.0
    es = math.exp(s)
    eu = math.exp(u)
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * math.exp(beta * u - c * es * (eu - 1.0))
    if sigma:
        total *= math.exp(0.5 * sigma * omega * (s + u))
    return total


def part_kernel(s, t, pcs, nodes, weights, which, omega=0.0, sigma=0):
    total = 0.0
    for u, w in zip(nodes, weights):
        left = a_value(s, float(u), pcs, omega, sigma)
        right = a_value(t, float(u), pcs, omega, sigma)
        if which == "h":
            factor = 1.0
        elif which == "u":
            factor = float(u)
        elif which == "s+t":
            factor = 0.5 * (s + t)
        elif which == "full-linear":
            factor = float(u) + 0.5 * (s + t)
        elif which == "full-cosh":
            center = float(u) + 0.5 * (s + t)
            factor = center * math.cosh(omega * center)
        elif which == "min-linear":
            factor = float(u) + min(s, t)
        elif which == "abs-linear":
            factor = 0.5 * abs(s - t)
        else:
            raise ValueError(which)
        total += float(w) * factor * left * right
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--vmax", type=float, default=8.0)
    parser.add_argument("--n", type=int, default=60)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--quad", type=int, default=220)
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--sigma", type=int, choices=(-1, 0, 1), default=0)
    args = parser.parse_args()

    pcs = pieces(args.kind)
    leg_nodes, leg_weights = np.polynomial.legendre.leggauss(args.quad)
    nodes = 0.5 * args.umax * (leg_nodes + 1.0)
    weights = 0.5 * args.umax * leg_weights
    pts = np.linspace(0.0, args.vmax, args.n)
    print(
        f"reduced parts kind={args.kind} omega={args.omega:g} sigma={args.sigma} "
        f"v=[0,{args.vmax:g}] n={args.n} u=[0,{args.umax:g}] quad={args.quad}"
    )
    for which in ("h", "u", "s+t", "min-linear", "abs-linear", "full-linear", "full-cosh"):
        mat = np.zeros((args.n, args.n), dtype=float)
        for i, s in enumerate(pts):
            for j, t in enumerate(pts[: i + 1]):
                value = part_kernel(float(s), float(t), pcs, nodes, weights, which, args.omega, args.sigma)
                mat[i, j] = mat[j, i] = value
        vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
        print(
            f"  {which}: min={vals[0]:.12e} max={vals[-1]:.12e} "
            f"neg={(vals < -1e-12).sum()}"
        )
        print("    low=", " ".join(f"{v:.4e}" for v in vals[:8]))


if __name__ == "__main__":
    main()
