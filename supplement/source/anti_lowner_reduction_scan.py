#!/usr/bin/env python3
import argparse
import math

from second_order_cosh_ibp import finite_modes

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "anti_lowner_reduction_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


PI = math.pi


def score_a(v: float, modes):
    logs = []
    ev = math.exp(v)
    for n, weight in modes:
        logs.append(math.log(weight) - PI * n * n * ev)
    shift = max(logs)
    z = 0.0
    mean_c = 0.0
    for n, weight in modes:
        c = PI * n * n
        p = math.exp(math.log(weight) - c * ev - shift)
        z += p
        mean_c += p * c
    return ev * mean_c / z - 0.25


def reduced_kernel(s: float, t: float, modes, alpha: float):
    # alpha corresponds to exp(alpha(s+t)); it is a diagonal congruence.
    if s == 0.0 and t == 0.0:
        return 0.0
    return (s + t) * math.exp(alpha * (s + t)) / (score_a(s, modes) + score_a(t, modes))


def matrix(points, modes, alpha: float, normalize: bool):
    n = len(points)
    mat = np.zeros((n, n), dtype=float)
    for i, s in enumerate(points):
        for j, t in enumerate(points[: i + 1]):
            value = reduced_kernel(float(s), float(t), modes, alpha)
            mat[i, j] = mat[j, i] = value
    if normalize:
        diag = np.sqrt(np.maximum(np.abs(np.diag(mat)), 1e-300))
        mat = mat / diag[:, None] / diag[None, :]
    return mat


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--alpha", type=float, default=0.0)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--tol", type=float, default=1e-10)
    args = parser.parse_args()

    modes = finite_modes(args.kind)
    points = points_for(args)
    mat = matrix(points, modes, args.alpha, args.normalize)
    vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
    print(
        f"kind={args.kind} grid={args.grid} v=[{args.vmin:g},{args.vmax:g}] "
        f"n={args.n} alpha={args.alpha:g} normalize={args.normalize}"
    )
    print(
        f"  inertia: neg={(vals < -args.tol).sum()} "
        f"zero={(abs(vals) <= args.tol).sum()} pos={(vals > args.tol).sum()}"
    )
    print(f"  min={vals[0]:.12e} max={vals[-1]:.12e}")
    print("  low=", " ".join(f"{v:.4e}" for v in vals[: min(12, len(vals))]))


if __name__ == "__main__":
    main()
