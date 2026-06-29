#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_symmetric_form.py requires numpy; run with python") from exc


def kappa(c, x, y):
    xy = x * y
    return xy * (c * xy - 1.5) * math.exp(-c * xy)


def eigvals(mat):
    return np.linalg.eigvalsh((mat + mat.T) / 2.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--Xmax", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=120)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    pts, wts = quadrature(args.Xmax - 1.0, args.order)
    pts = pts + 1.0
    mu_wts = wts * pts ** -0.5
    root = np.sqrt(mu_wts)
    kmat = np.zeros((args.order, args.order), dtype=float)
    for i, x in enumerate(pts):
        for j, y in enumerate(pts[: i + 1]):
            value = root[i] * kappa(args.c, float(x), float(y)) * root[j]
            kmat[i, j] = kmat[j, i] = value

    logm = np.diag(np.log(pts))
    # Quadratic form <Kf,LKf> + <Kf,K L f>.
    target = kmat @ logm @ kmat + 0.5 * (kmat @ kmat @ logm + logm @ kmat @ kmat)

    kvals = eigvals(kmat)
    tvals = eigvals(target)
    print(
        f"endpoint symmetric form X=[1,{args.Xmax:g}] order={args.order} c={args.c:g}"
    )
    print(
        f"  K min={kvals[0]:.12e} max={kvals[-1]:.12e} "
        f"neg={(kvals < -args.tol).sum()}"
    )
    print(
        f"  anti-commutator min={tvals[0]:.12e} max={tvals[-1]:.12e} "
        f"neg={(tvals < -args.tol).sum()}"
    )
    print("  low K eigenvalues:", " ".join(f"{v:.3e}" for v in kvals[:8]))
    print("  low target eigenvalues:", " ".join(f"{v:.3e}" for v in tvals[:8]))


if __name__ == "__main__":
    main()
