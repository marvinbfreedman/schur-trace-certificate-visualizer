#!/usr/bin/env python3
import argparse
import math

from anti_lowner_reduction_scan import score_a
from second_order_cosh_ibp import finite_modes

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "boundary_mode_fit.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def points_and_weights(args):
    if args.grid == "linear":
        pts = np.linspace(args.vmin, args.vmax, args.n)
    elif args.grid == "quadratic":
        z = np.linspace(0.0, 1.0, args.n)
        pts = args.vmin + (args.vmax - args.vmin) * z * z
    elif args.grid == "geometric":
        if args.vmin <= 0:
            raise ValueError("geometric grid requires --vmin > 0")
        pts = np.geomspace(args.vmin, args.vmax, args.n)
    else:
        raise ValueError(args.grid)
    weights = np.zeros_like(pts)
    weights[0] = 0.5 * (pts[1] - pts[0])
    weights[-1] = 0.5 * (pts[-1] - pts[-2])
    weights[1:-1] = 0.5 * (pts[2:] - pts[:-2])
    return pts, weights


def boundary_value(s, t, omega, modes):
    if s == 0.0 and t == 0.0:
        return 0.0
    u = 0.5 * (s + t)
    return u * math.cosh(omega * u) / (score_a(s, modes) + score_a(t, modes))


def weighted_matrix(pts, weights, omega, modes):
    n = len(pts)
    mat = np.zeros((n, n), dtype=float)
    root_w = np.sqrt(weights)
    for i, s in enumerate(pts):
        for j, t in enumerate(pts[: i + 1]):
            value = root_w[i] * boundary_value(float(s), float(t), omega, modes) * root_w[j]
            mat[i, j] = mat[j, i] = value
    return mat


def normalize_weighted(vec, weights):
    norm = math.sqrt(float(np.sum(weights * vec * vec)))
    return vec / norm if norm else vec


def weighted_basis(pts, weights):
    raw = [
        ("1", np.ones_like(pts)),
        ("s", pts.copy()),
        ("s^2", pts * pts),
        ("s^3", pts * pts * pts),
        ("exp(-s)", np.exp(-pts)),
        ("s exp(-s)", pts * np.exp(-pts)),
        ("endpoint", np.exp(-20.0 * pts)),
    ]
    out = []
    for name, vec in raw:
        candidate = vec.copy()
        for _, prev in out:
            candidate -= np.sum(weights * candidate * prev) * prev
        norm = math.sqrt(float(np.sum(weights * candidate * candidate)))
        if norm > 1e-14:
            out.append((name, candidate / norm))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="raw1")
    parser.add_argument("--omega", type=float, default=0.0)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=160)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--modes", type=int, default=5)
    args = parser.parse_args()

    modes = finite_modes(args.kind)
    pts, weights = points_and_weights(args)
    mat = weighted_matrix(pts, weights, args.omega, modes)
    vals, vecs = np.linalg.eigh((mat + mat.T) / 2.0)
    basis = weighted_basis(pts, weights)
    print(
        f"kind={args.kind} omega={args.omega:g} grid={args.grid} "
        f"v=[{args.vmin:g},{args.vmax:g}] n={args.n}"
    )
    for k in range(min(args.modes, len(vals))):
        eig = vals[k]
        fn = vecs[:, k] / np.sqrt(weights)
        fn = normalize_weighted(fn, weights)
        print(f"mode {k + 1}: eig={eig:.12e}")
        corr = []
        for name, bvec in basis:
            corr.append((abs(float(np.sum(weights * fn * bvec))), name, float(np.sum(weights * fn * bvec))))
        corr.sort(reverse=True)
        print("  correlations:", " ".join(f"{name}={signed:+.4f}" for _, name, signed in corr[:5]))
        sample_idx = np.linspace(0, len(pts) - 1, 8, dtype=int)
        print(
            "  samples:",
            " ".join(f"({pts[i]:.3g},{fn[i]:+.3e})" for i in sample_idx),
        )


if __name__ == "__main__":
    main()
