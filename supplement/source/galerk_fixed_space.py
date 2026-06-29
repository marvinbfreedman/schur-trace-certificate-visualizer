#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import boundary_red, full_red, pieces
from second_order_cosh_ibp import finite_modes

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "galerk_fixed_space.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def quadrature(interval_end, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    pts = 0.5 * interval_end * (nodes + 1.0)
    wts = 0.5 * interval_end * weights
    return pts, wts


def weighted_matrices(kind, omega, interval_end, quad_order, laguerre_order):
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
    return pts, wts, boundary, full, full - boundary


def orthonormalize(cols):
    mat = np.column_stack(cols)
    q, r = np.linalg.qr(mat)
    keep = np.abs(np.diag(r)) > 1e-12
    return q[:, keep]


def trial_space(pts, wts, interval_end, name):
    root = np.sqrt(wts)
    if name == "center_exp5":
        funcs = [pts - 0.5 * interval_end, np.exp(-5.0 * pts)]
    elif name == "s_exp1":
        funcs = [pts, np.exp(-pts)]
    elif name == "one_exp5":
        funcs = [np.ones_like(pts), np.exp(-5.0 * pts)]
    elif name == "center_exp1":
        funcs = [pts - 0.5 * interval_end, np.exp(-pts)]
    else:
        raise ValueError(name)
    return orthonormalize([root * fn for fn in funcs])


def complete_basis(q, n):
    full, _ = np.linalg.qr(q, mode="complete")
    return full[:, q.shape[1] : n]


def generalized_ratio(pos, neg_basis, neg_vals):
    if neg_basis.shape[1] == 0:
        return None
    block = neg_basis.T @ ((pos + pos.T) / 2.0) @ neg_basis
    scale = np.diag(1.0 / np.sqrt(-neg_vals))
    return np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)


def analyze(boundary, tail, full, q, tol):
    n = full.shape[0]
    rest = complete_basis(q, n)
    kk = (full + full.T) / 2.0
    aa = q.T @ kk @ q
    ab = q.T @ kk @ rest
    bb = rest.T @ kk @ rest
    bb_vals, bb_vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = bb_vals > tol
    inv_bb = (bb_vecs[:, keep] * (1.0 / bb_vals[keep])) @ bb_vecs[:, keep].T
    schur = aa - ab @ inv_bb @ ab.T
    range_resid = np.linalg.norm(ab @ (np.eye(len(bb_vals)) - bb_vecs[:, keep] @ bb_vecs[:, keep].T))

    cb = rest.T @ ((boundary + boundary.T) / 2.0) @ rest
    tb = rest.T @ ((tail + tail.T) / 2.0) @ rest
    c_vals, c_vecs = np.linalg.eigh((cb + cb.T) / 2.0)
    neg = c_vals < -tol
    ratios = generalized_ratio(tb, c_vecs[:, neg], c_vals[neg])
    return {
        "full_min": np.linalg.eigvalsh(kk)[0],
        "tail_min": np.linalg.eigvalsh((tail + tail.T) / 2.0)[0],
        "bb_min": bb_vals[0],
        "bb_pos_rank": int(keep.sum()),
        "bb_min_pos": float(bb_vals[keep][0]) if keep.any() else 0.0,
        "range_resid": range_resid,
        "schur_min": np.linalg.eigvalsh((schur + schur.T) / 2.0)[0],
        "res_neg_count": int(neg.sum()),
        "res_neg_min": c_vals[0],
        "ratio_min": None if ratios is None else ratios[0],
        "ratio_max": None if ratios is None else ratios[-1],
        "boundary_lows": np.linalg.eigvalsh((boundary + boundary.T) / 2.0)[:8],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--quad", type=int, nargs="+", default=[40, 60, 80])
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument(
        "--space",
        choices=("center_exp5", "s_exp1", "one_exp5", "center_exp1"),
        default="center_exp5",
    )
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    print(
        f"galerkin fixed-space kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} laguerre={args.laguerre} space={args.space}"
    )
    for order in args.quad:
        pts, wts, boundary, full, tail = weighted_matrices(
            args.kind, args.omega, args.L, order, args.laguerre
        )
        q = trial_space(pts, wts, args.L, args.space)
        data = analyze(boundary, tail, full, q, args.tol)
        ratio_text = "none" if data["ratio_min"] is None else f"{data['ratio_min']:.12e}"
        print(f"quad={order}")
        print("  boundary low=", " ".join(f"{v:.4e}" for v in data["boundary_lows"]))
        print(
            f"  full_min={data['full_min']:.12e} tail_min={data['tail_min']:.12e} "
            f"bb_min={data['bb_min']:.12e}"
        )
        print(
            f"  schur_min={data['schur_min']:.12e} "
            f"bb_pos_rank={data['bb_pos_rank']} bb_min_pos={data['bb_min_pos']:.12e} "
            f"range_resid={data['range_resid']:.12e}"
        )
        print(
            f"  residual: neg={data['res_neg_count']} "
            f"min={data['res_neg_min']:.12e} tail/res={ratio_text}"
        )


if __name__ == "__main__":
    main()
