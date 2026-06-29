#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import boundary_red, full_red, pieces
from second_order_cosh_ibp import finite_modes

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "legendre_certificate.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def shifted_legendre_values(x, order, interval_end):
    z = 2.0 * x / interval_end - 1.0
    values = np.zeros((len(x), order), dtype=float)
    values[:, 0] = 1.0 / math.sqrt(interval_end)
    if order == 1:
        return values
    p0 = np.ones_like(z)
    p1 = z.copy()
    values[:, 1] = math.sqrt(3.0 / interval_end) * p1
    for n in range(1, order - 1):
        p2 = ((2 * n + 1) * z * p1 - n * p0) / (n + 1)
        values[:, n + 1] = math.sqrt((2 * n + 3) / interval_end) * p2
        p0, p1 = p1, p2
    return values


def quadrature(interval_end, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    pts = 0.5 * interval_end * (nodes + 1.0)
    wts = 0.5 * interval_end * weights
    return pts, wts


def weighted_kernel(kind, omega, interval_end, quad_order, laguerre_order):
    modes = finite_modes(kind)
    pcs = pieces(kind)
    pts, wts = quadrature(interval_end, quad_order)
    lag_nodes, lag_weights = np.polynomial.laguerre.laggauss(laguerre_order)
    root = np.sqrt(wts)
    boundary = np.zeros((quad_order, quad_order), dtype=float)
    full = np.zeros((quad_order, quad_order), dtype=float)
    for i, s in enumerate(pts):
        for j, t in enumerate(pts[: i + 1]):
            b = boundary_red(float(s), float(t), omega, modes)
            k = full_red(float(s), float(t), omega, pcs, lag_nodes, lag_weights)
            boundary[i, j] = boundary[j, i] = root[i] * b * root[j]
            full[i, j] = full[j, i] = root[i] * k * root[j]
    return pts, wts, boundary, full


def project_to_legendre(weighted_mat, pts, wts, basis_order, interval_end):
    vals = shifted_legendre_values(pts, basis_order, interval_end)
    weighted_vals = np.sqrt(wts)[:, None] * vals
    return weighted_vals.T @ weighted_mat @ weighted_vals


def trial_coefficients(pts, wts, basis_order, interval_end, trial):
    vals = shifted_legendre_values(pts, basis_order, interval_end)
    if trial == "center_exp5":
        funcs = [pts - 0.5 * interval_end, np.exp(-5.0 * pts)]
    elif trial == "s_exp1":
        funcs = [pts, np.exp(-pts)]
    elif trial == "one_exp5":
        funcs = [np.ones_like(pts), np.exp(-5.0 * pts)]
    else:
        raise ValueError(trial)
    coeffs = []
    for fn in funcs:
        coeffs.append(vals.T @ (wts * fn))
    mat = np.column_stack(coeffs)
    q, r = np.linalg.qr(mat)
    keep = np.abs(np.diag(r)) > 1e-12
    return q[:, keep]


def complement(q, size):
    full, _ = np.linalg.qr(q, mode="complete")
    return full[:, q.shape[1] : size]


def schur_data(kmat, emat, tol):
    rest = complement(emat, kmat.shape[0])
    kk = (kmat + kmat.T) / 2.0
    aa = emat.T @ kk @ emat
    ab = emat.T @ kk @ rest
    bb = rest.T @ kk @ rest
    vals, vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = vals > tol
    inv_bb = (vecs[:, keep] * (1.0 / vals[keep])) @ vecs[:, keep].T
    sc = aa - ab @ inv_bb @ ab.T
    resid = np.linalg.norm(ab @ (np.eye(len(vals)) - vecs[:, keep] @ vecs[:, keep].T))
    return rest, vals, np.linalg.eigvalsh((sc + sc.T) / 2.0), resid


def residual_ratio(cmat, tmat, rest, tol):
    cb = rest.T @ ((cmat + cmat.T) / 2.0) @ rest
    tb = rest.T @ ((tmat + tmat.T) / 2.0) @ rest
    c_vals, c_vecs = np.linalg.eigh((cb + cb.T) / 2.0)
    neg = c_vals < -tol
    if not neg.any():
        return c_vals, None
    block = c_vecs[:, neg].T @ tb @ c_vecs[:, neg]
    scale = np.diag(1.0 / np.sqrt(-c_vals[neg]))
    ratios = np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)
    return c_vals, ratios


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--basis", type=int, nargs="+", default=[8, 12, 16, 20])
    parser.add_argument("--quad", type=int, default=160)
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument("--trial", choices=("center_exp5", "s_exp1", "one_exp5"), default="center_exp5")
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    pts, wts, boundary_w, full_w = weighted_kernel(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    print(
        f"legendre certificate kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"quad={args.quad} laguerre={args.laguerre} trial={args.trial}"
    )
    for order in args.basis:
        cmat = project_to_legendre(boundary_w, pts, wts, order, args.L)
        kmat = project_to_legendre(full_w, pts, wts, order, args.L)
        tmat = kmat - cmat
        emat = trial_coefficients(pts, wts, order, args.L, args.trial)
        rest, bb_vals, sc_vals, range_resid = schur_data(kmat, emat, args.tol)
        c_vals, ratios = residual_ratio(cmat, tmat, rest, args.tol)
        kvals = np.linalg.eigvalsh((kmat + kmat.T) / 2.0)
        print(f"basis={order}")
        print(
            f"  K_min={kvals[0]:.12e} K_rest_min={bb_vals[0]:.12e} "
            f"Schur_min={sc_vals[0]:.12e} range_resid={range_resid:.12e}"
        )
        print(
            f"  C_rest_neg={(c_vals < -args.tol).sum()} "
            f"C_rest_min={c_vals[0]:.12e}"
        )
        if ratios is not None:
            print(
                f"  tail/residual ratio: min={ratios[0]:.12e} "
                f"max={ratios[-1]:.12e}"
            )


if __name__ == "__main__":
    main()
