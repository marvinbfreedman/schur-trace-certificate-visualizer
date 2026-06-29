#!/usr/bin/env python3
import argparse

from reduced_exact_finite import matrix

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "split_scale_schur.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def orthogonal_completion(basis, size):
    if basis.shape[1] == 0:
        return np.eye(size)
    q, _ = np.linalg.qr(basis, mode="complete")
    return q


def schur(kernel, lead_basis, tol):
    n = kernel.shape[0]
    q = orthogonal_completion(lead_basis, n)
    k = lead_basis.shape[1]
    left = q[:, :k]
    rest = q[:, k:]
    kk = (kernel + kernel.T) / 2.0
    aa = left.T @ kk @ left
    ab = left.T @ kk @ rest
    bb = rest.T @ kk @ rest
    vals, vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = vals > tol
    inv_bb = (vecs[:, keep] * (1.0 / vals[keep])) @ vecs[:, keep].T
    complement = aa - ab @ inv_bb @ ab.T
    range_residual = np.linalg.norm(ab @ (np.eye(len(vals)) - vecs[:, keep] @ vecs[:, keep].T))
    return vals, np.linalg.eigvalsh((complement + complement.T) / 2.0), range_residual


def generalized_ratio(pos, neg_basis, neg_vals):
    if neg_basis.shape[1] == 0:
        return np.array([])
    block = neg_basis.T @ ((pos + pos.T) / 2.0) @ neg_basis
    scale = np.diag(1.0 / np.sqrt(-neg_vals))
    return np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--laguerre", type=int, default=120)
    parser.add_argument("--lead", type=int, default=2)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    boundary = matrix(args, "boundary")
    tail = matrix(args, "tail")
    full = boundary + tail
    bvals, bvecs = np.linalg.eigh((boundary + boundary.T) / 2.0)
    lead_basis = bvecs[:, : args.lead]
    rest_boundary = bvecs[:, args.lead :]
    rest_bvals = bvals[args.lead :]
    residual_neg_mask = rest_bvals < -args.tol
    residual_neg_basis = rest_boundary[:, residual_neg_mask]
    residual_neg_vals = rest_bvals[residual_neg_mask]

    bb_vals, sc_vals, range_resid = schur(full, lead_basis, args.tol)
    ratios = generalized_ratio(tail, residual_neg_basis, residual_neg_vals)
    fvals = np.linalg.eigvalsh((full + full.T) / 2.0)
    tvals = np.linalg.eigvalsh((tail + tail.T) / 2.0)

    print(
        f"split-scale kind={args.kind} omega={args.omega:g} grid={args.grid} "
        f"v=[{args.vmin:g},{args.vmax:g}] n={args.n} lead={args.lead}"
    )
    print("  boundary low=", " ".join(f"{v:.6e}" for v in bvals[: min(8, len(bvals))]))
    print(f"  full min={fvals[0]:.12e} tail min={tvals[0]:.12e}")
    print(
        f"  lead Schur: rest_pos_rank={(bb_vals > args.tol).sum()} "
        f"rest_min_pos={bb_vals[bb_vals > args.tol][0]:.12e} "
        f"range_resid={range_resid:.12e} schur_min={sc_vals[0]:.12e}"
    )
    print(
        f"  residual negatives after lead: count={len(residual_neg_vals)} "
        f"low={' '.join(f'{v:.4e}' for v in residual_neg_vals[:6])}"
    )
    if len(ratios):
        print(
            f"  tail/residual(-C): min={ratios[0]:.12e} "
            f"max={ratios[-1]:.12e}"
        )


if __name__ == "__main__":
    main()
