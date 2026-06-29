#!/usr/bin/env python3
import argparse
import math

from second_order_cosh_ibp import finite_modes, scalar_base_values

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "reduced_exact_finite.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


PI = math.pi


def pieces(kind):
    out = []
    for n, weight in finite_modes(kind):
        c = PI * n * n
        out.append((weight * 4.0 * c * c, 2.25, c))
        out.append((weight * -6.0 * c, 1.25, c))
    return out


def psi(v, pcs):
    ev = math.exp(v)
    total = 0.0
    for coeff, beta, c in pcs:
        exponent = beta * v - c * ev
        if exponent > -745.0:
            total += coeff * math.exp(exponent)
    return total


def endpoint_ratios(v, pcs):
    denom = psi(v, pcs)
    ev = math.exp(v)
    ratios = []
    for coeff, beta, c in pcs:
        exponent = beta * v - c * ev
        value = coeff * math.exp(exponent) / denom if exponent > -745.0 else 0.0
        ratios.append((value, beta, c))
    return ratios


def boundary_red(s, t, omega, modes):
    center = 0.5 * (s + t)
    if center == 0.0:
        return 0.0
    _, a_s, _, _, _ = scalar_base_values(s, modes)
    _, a_t, _, _, _ = scalar_base_values(t, modes)
    return center * math.cosh(omega * center) / (a_s + a_t)


def laguerre_integral(alpha, p, center, omega, nodes, weights):
    # exp(alpha) * int_1^inf (center+log q) cosh(omega(center+log q))
    #   * q^p exp(-alpha q) dq
    # = 1/alpha * int_0^inf (...) exp(-r) dr, q=1+r/alpha.
    q = 1.0 + nodes / alpha
    u = center + np.log(q)
    return float(np.sum(weights * u * np.cosh(omega * u) * np.power(q, p)) / alpha)


def full_red(s, t, omega, pcs, nodes, weights):
    center = 0.5 * (s + t)
    left = endpoint_ratios(s, pcs)
    right = endpoint_ratios(t, pcs)
    es = math.exp(s)
    et = math.exp(t)
    total = 0.0
    for ratio_i, beta_i, c_i in left:
        if ratio_i == 0.0:
            continue
        for ratio_j, beta_j, c_j in right:
            if ratio_j == 0.0:
                continue
            alpha = c_i * es + c_j * et
            p = beta_i + beta_j - 1.0
            total += ratio_i * ratio_j * laguerre_integral(
                alpha, p, center, omega, nodes, weights
            )
    return total


def points_for(args):
    if args.grid == "linear":
        return np.linspace(args.vmin, args.vmax, args.n)
    if args.grid == "quadratic":
        z = np.linspace(0.0, 1.0, args.n)
        return args.vmin + (args.vmax - args.vmin) * z * z
    if args.grid == "geometric":
        if args.vmin <= 0:
            raise ValueError("geometric grid requires --vmin > 0")
        return np.geomspace(args.vmin, args.vmax, args.n)
    raise ValueError(args.grid)


def matrix(args, which):
    modes = finite_modes(args.kind)
    pcs = pieces(args.kind)
    nodes, weights = np.polynomial.laguerre.laggauss(args.laguerre)
    pts = points_for(args)
    out = np.zeros((args.n, args.n), dtype=float)
    for i, s in enumerate(pts):
        for j, t in enumerate(pts[: i + 1]):
            b = boundary_red(float(s), float(t), args.omega, modes)
            if which == "boundary":
                value = b
            elif which == "full":
                value = full_red(float(s), float(t), args.omega, pcs, nodes, weights)
            elif which == "tail":
                value = full_red(float(s), float(t), args.omega, pcs, nodes, weights) - b
            else:
                raise ValueError(which)
            out[i, j] = out[j, i] = value
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=60)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--laguerre", type=int, default=120)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    print(
        f"reduced exact finite kind={args.kind} omega={args.omega:g} "
        f"grid={args.grid} v=[{args.vmin:g},{args.vmax:g}] n={args.n} "
        f"laguerre={args.laguerre}"
    )
    mats = {name: matrix(args, name) for name in ("boundary", "tail", "full")}
    for name, mat in mats.items():
        vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
        print(
            f"  {name}: neg={(vals < -args.tol).sum()} "
            f"zero={(abs(vals) <= args.tol).sum()} pos={(vals > args.tol).sum()} "
            f"min={vals[0]:.12e} max={vals[-1]:.12e}"
        )
        print("    low=", " ".join(f"{v:.4e}" for v in vals[: min(10, len(vals))]))


if __name__ == "__main__":
    main()
