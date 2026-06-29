#!/usr/bin/env python3
import argparse
import math

from galerk_fixed_space import complete_basis, generalized_ratio, weighted_matrices

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("optimize_fixed_space.py requires numpy; run with python, not python3") from exc


def orthonormalize(cols):
    mat = np.column_stack(cols)
    q, r = np.linalg.qr(mat)
    keep = np.abs(np.diag(r)) > 1e-12
    return q[:, keep]


def trial_space(pts, wts, interval_end, beta, affine):
    root = np.sqrt(wts)
    if affine == "center":
        first = pts - 0.5 * interval_end
    elif affine == "one":
        first = np.ones_like(pts)
    elif affine == "s":
        first = pts.copy()
    else:
        raise ValueError(affine)
    return orthonormalize([root * first, root * np.exp(-beta * pts)])


def schur_stats(boundary, tail, full, q, tol):
    n = full.shape[0]
    rest = complete_basis(q, n)
    kk = (full + full.T) / 2.0
    aa = q.T @ kk @ q
    ab = q.T @ kk @ rest
    bb = rest.T @ kk @ rest
    bb_vals, bb_vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = bb_vals > tol
    if keep.any():
        inv_bb = (bb_vecs[:, keep] * (1.0 / bb_vals[keep])) @ bb_vecs[:, keep].T
        sc = aa - ab @ inv_bb @ ab.T
        projector = bb_vecs[:, keep] @ bb_vecs[:, keep].T
        range_resid = np.linalg.norm(ab @ (np.eye(len(bb_vals)) - projector))
        rest_min_pos = float(bb_vals[keep][0])
    else:
        sc = aa
        range_resid = float(np.linalg.norm(ab))
        rest_min_pos = 0.0

    cb = rest.T @ ((boundary + boundary.T) / 2.0) @ rest
    tb = rest.T @ ((tail + tail.T) / 2.0) @ rest
    c_vals, c_vecs = np.linalg.eigh((cb + cb.T) / 2.0)
    neg = c_vals < -tol
    ratios = generalized_ratio(tb, c_vecs[:, neg], c_vals[neg])

    b_vals, b_vecs = np.linalg.eigh((boundary + boundary.T) / 2.0)
    lead = b_vecs[:, :2]
    capture = np.linalg.svd(q.T @ lead, compute_uv=False)

    return {
        "schur_min": float(np.linalg.eigvalsh((sc + sc.T) / 2.0)[0]),
        "rest_min": float(bb_vals[0]),
        "rest_min_pos": rest_min_pos,
        "rest_pos_rank": int(keep.sum()),
        "range_resid": float(range_resid),
        "res_neg_count": int(neg.sum()),
        "res_neg_min": float(c_vals[0]),
        "ratio_min": None if ratios is None else float(ratios[0]),
        "ratio_max": None if ratios is None else float(ratios[-1]),
        "capture_min": float(capture[-1]),
        "capture_max": float(capture[0]),
        "boundary_low": b_vals[:4],
    }


def beta_values(args):
    if args.betas:
        return args.betas
    return np.linspace(args.beta_min, args.beta_max, args.beta_n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--quad", type=int, default=120)
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument("--affine", choices=("center", "one", "s"), default="center")
    parser.add_argument("--beta-min", type=float, default=1.0)
    parser.add_argument("--beta-max", type=float, default=12.0)
    parser.add_argument("--beta-n", type=int, default=45)
    parser.add_argument("--betas", type=float, nargs="*", default=None)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    pts, wts, boundary, full, tail = weighted_matrices(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    print(
        f"optimize fixed E kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"quad={args.quad} laguerre={args.laguerre} affine={args.affine}"
    )
    rows = []
    for beta in beta_values(args):
        q = trial_space(pts, wts, args.L, float(beta), args.affine)
        if q.shape[1] != 2:
            continue
        data = schur_stats(boundary, tail, full, q, args.tol)
        ratio = data["ratio_min"] if data["ratio_min"] is not None else float("inf")
        score = min(data["schur_min"] / 1e-7, ratio)
        rows.append((score, float(beta), data))

    rows.sort(key=lambda row: (-row[0], -row[2]["schur_min"], -row[2]["capture_min"]))
    for _, beta, data in rows[: args.top]:
        ratio = "none" if data["ratio_min"] is None else f"{data['ratio_min']:.6e}"
        print(
            f"beta={beta:.6g} schur={data['schur_min']:.6e} "
            f"Krest={data['rest_min']:.3e} res_neg={data['res_neg_count']} "
            f"Cres_min={data['res_neg_min']:.6e} tail/res={ratio} "
            f"capture=({data['capture_min']:.6f},{data['capture_max']:.6f}) "
            f"range_resid={data['range_resid']:.3e}"
        )
    if rows:
        best = rows[0][2]
        lows = " ".join(f"{x:.4e}" for x in best["boundary_low"])
        print(f"boundary lows for best beta: {lows}")


if __name__ == "__main__":
    main()
