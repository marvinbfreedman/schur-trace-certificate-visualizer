#!/usr/bin/env python3
import argparse
import itertools
import math

from reduced_exact_finite import matrix, points_for

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "fixed_space_schur.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def candidates(pts, length):
    x = pts / length
    return {
        "1": np.ones_like(pts),
        "s": pts,
        "s2": pts * pts,
        "s3": pts * pts * pts,
        "center": pts - 0.5 * length,
        "quad": (pts - 0.5 * length) ** 2,
        "cubic": (pts - 0.5 * length) ** 3,
        "exp1": np.exp(-pts),
        "sexp1": pts * np.exp(-pts),
        "exp5": np.exp(-5.0 * pts),
        "exp10": np.exp(-10.0 * pts),
        "exp20": np.exp(-20.0 * pts),
        "x1x": x * (1.0 - x),
        "x1x2": x * (1.0 - x) * (1.0 - 2.0 * x),
    }


def orthonormalize(cols):
    mat = np.column_stack(cols)
    q, r = np.linalg.qr(mat)
    keep = np.abs(np.diag(r)) > 1e-12
    return q[:, keep]


def complement_basis(q, n):
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
    rest = complement_basis(q, n)
    kk = (full + full.T) / 2.0
    aa = q.T @ kk @ q
    ab = q.T @ kk @ rest
    bb = rest.T @ kk @ rest
    bb_vals, bb_vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = bb_vals > tol
    if keep.any():
        inv_bb = (bb_vecs[:, keep] * (1.0 / bb_vals[keep])) @ bb_vecs[:, keep].T
        schur_vals = np.linalg.eigvalsh((aa - ab @ inv_bb @ ab.T + (aa - ab @ inv_bb @ ab.T).T) / 2.0)
        resid = np.linalg.norm(ab @ (np.eye(len(bb_vals)) - bb_vecs[:, keep] @ bb_vecs[:, keep].T))
        min_pos = bb_vals[keep][0]
    else:
        schur_vals = np.linalg.eigvalsh((aa + aa.T) / 2.0)
        resid = np.linalg.norm(ab)
        min_pos = 0.0

    cb = rest.T @ ((boundary + boundary.T) / 2.0) @ rest
    c_vals, c_vecs = np.linalg.eigh((cb + cb.T) / 2.0)
    neg = c_vals < -tol
    neg_basis = rest @ c_vecs[:, neg]
    ratios = generalized_ratio(tail, neg_basis, c_vals[neg])
    return {
        "rest_min": bb_vals[0],
        "rest_min_pos": min_pos,
        "rest_pos_rank": int(keep.sum()),
        "range_resid": resid,
        "schur_min": schur_vals[0],
        "res_neg_count": int(neg.sum()),
        "res_neg_min": c_vals[0],
        "ratio_min": None if ratios is None else ratios[0],
        "ratio_max": None if ratios is None else ratios[-1],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--laguerre", type=int, default=120)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--pair", default="")
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    pts = points_for(args)
    args.n = len(pts)
    boundary = matrix(args, "boundary")
    tail = matrix(args, "tail")
    full = boundary + tail
    cand = candidates(pts, args.vmax - args.vmin)

    if args.pair:
        pairs = [tuple(args.pair.split(","))]
    else:
        names = list(cand)
        pairs = list(itertools.combinations(names, 2))

    rows = []
    for pair in pairs:
        if pair[0] not in cand or pair[1] not in cand:
            raise SystemExit(f"unknown pair {pair}")
        q = orthonormalize([cand[pair[0]], cand[pair[1]]])
        if q.shape[1] != 2:
            continue
        data = analyze(boundary, tail, full, q, args.tol)
        rows.append((data["schur_min"], data["ratio_min"] if data["ratio_min"] is not None else float("inf"), pair, data))

    rows.sort(key=lambda row: (-row[0], -(row[1] if row[1] is not None else -1.0)))
    print(
        f"fixed-space kind={args.kind} omega={args.omega:g} grid={args.grid} "
        f"v=[{args.vmin:g},{args.vmax:g}] n={args.n}"
    )
    for schur_min, ratio_min, pair, data in rows[: args.top]:
        ratio_text = "none" if data["ratio_min"] is None else f"{data['ratio_min']:.6e}"
        print(
            f"  {pair[0]},{pair[1]}: schur={data['schur_min']:.6e} "
            f"rest_min={data['rest_min']:.6e} rest_pos={data['rest_pos_rank']} "
            f"res_neg={data['res_neg_count']} res_min={data['res_neg_min']:.6e} "
            f"tail/res={ratio_text} range_resid={data['range_resid']:.3e}"
        )


if __name__ == "__main__":
    main()
