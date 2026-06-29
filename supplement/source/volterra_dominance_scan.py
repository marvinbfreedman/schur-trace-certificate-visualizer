#!/usr/bin/env python3
import argparse

from second_order_cosh_ibp import (
    direct_integrated_matrix,
    finite_modes,
    scalar_base_integrated_matrix,
    scalar_base_kernel_matrix,
)

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "volterra_dominance_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def split_symmetric(mat, tol):
    vals, vecs = np.linalg.eigh((mat + mat.T) / 2.0)
    neg_mask = vals < -tol
    pos_mask = vals > tol
    neg = (vecs[:, neg_mask] * (-vals[neg_mask])) @ vecs[:, neg_mask].T
    pos = (vecs[:, pos_mask] * vals[pos_mask]) @ vecs[:, pos_mask].T
    return vals, pos, neg, vecs[:, neg_mask], vals[neg_mask]


def restricted_eigs(mat, basis):
    if basis.shape[1] == 0:
        return np.array([])
    return np.linalg.eigvalsh((basis.T @ ((mat + mat.T) / 2.0) @ basis))


def generalized_tail_ratios(tail, neg_basis, neg_vals):
    if neg_basis.shape[1] == 0:
        return np.array([])
    # In the C-negative eigenspace, -C is diagonal with entries -neg_vals.
    t = neg_basis.T @ ((tail + tail.T) / 2.0) @ neg_basis
    scale = np.diag(1.0 / np.sqrt(-neg_vals))
    return np.linalg.eigvalsh((scale @ t @ scale + scale @ t.T @ scale) / 2.0)


def schur_data(kernel, neg_basis, tol):
    if neg_basis.shape[1] == 0:
        return None
    size = kernel.shape[0]
    q, _ = np.linalg.qr(neg_basis, mode="complete")
    rank_neg = neg_basis.shape[1]
    left = q[:, :rank_neg]
    rest = q[:, rank_neg:size]
    kk = (kernel + kernel.T) / 2.0
    nn = left.T @ kk @ left
    nr = left.T @ kk @ rest
    rr = rest.T @ kk @ rest
    vals, vecs = np.linalg.eigh((rr + rr.T) / 2.0)
    keep = vals > tol
    if not keep.any():
        return {
            "rr_rank": 0,
            "rr_min": 0.0,
            "range_residual": float(np.linalg.norm(nr)),
            "schur": np.linalg.eigvalsh((nn + nn.T) / 2.0),
        }
    projector = vecs[:, keep] @ vecs[:, keep].T
    inv_rr = (vecs[:, keep] * (1.0 / vals[keep])) @ vecs[:, keep].T
    schur = nn - nr @ inv_rr @ nr.T
    return {
        "rr_rank": int(keep.sum()),
        "rr_min": float(vals[keep][0]),
        "range_residual": float(np.linalg.norm(nr @ (np.eye(len(vals)) - projector))),
        "schur": np.linalg.eigvalsh((schur + schur.T) / 2.0),
    }


def scan(args):
    modes = finite_modes(args.kind)
    print(
        f"kind={args.kind} omega={args.omega:g} x=[{args.xmin:g},{args.xmax:g}] "
        f"u=[0,{args.umax:g}] intervals={args.intervals}"
    )
    for n in args.sizes:
        xs = np.linspace(args.xmin, args.xmax, n)
        boundary = scalar_base_kernel_matrix(xs, args.omega, modes, "C", 0.0)
        tail = scalar_base_integrated_matrix(
            xs, args.omega, modes, "D", args.umax, args.intervals
        )
        direct = direct_integrated_matrix(xs, args.omega, modes, args.umax, args.intervals)
        c_vals, c_pos, c_neg, neg_basis, neg_vals = split_symmetric(boundary, args.tol)
        tail_vals = np.linalg.eigvalsh((tail + tail.T) / 2.0)
        direct_vals = np.linalg.eigvalsh((direct + direct.T) / 2.0)
        residual_vals = np.linalg.eigvalsh((tail - c_neg + (tail - c_neg).T) / 2.0)
        ratios = generalized_tail_ratios(tail, neg_basis, neg_vals)
        tail_on_neg = restricted_eigs(tail, neg_basis)
        schur = schur_data(direct, neg_basis, args.tol)
        print(f"n={n}")
        print(
            f"  C inertia: neg={(c_vals < -args.tol).sum()} "
            f"zero={(abs(c_vals) <= args.tol).sum()} pos={(c_vals > args.tol).sum()} "
            f"min={c_vals[0]:.12e} max={c_vals[-1]:.12e}"
        )
        if len(neg_vals):
            print(
                f"  -C on C-neg: min={(-neg_vals).min():.12e} "
                f"max={(-neg_vals).max():.12e}"
            )
            print(
                f"  tail on C-neg: min={tail_on_neg[0]:.12e} "
                f"max={tail_on_neg[-1]:.12e}"
            )
            print(
                f"  generalized tail/(-C): min={ratios[0]:.12e} "
                f"max={ratios[-1]:.12e}"
            )
            print(
                f"  K Schur wrt C-neg: rr_rank={schur['rr_rank']} "
                f"rr_min={schur['rr_min']:.12e} "
                f"range_resid={schur['range_residual']:.12e} "
                f"schur_min={schur['schur'][0]:.12e}"
            )
        print(
            f"  tail min={tail_vals[0]:.12e} direct min={direct_vals[0]:.12e} "
            f"tail-Cneg min={residual_vals[0]:.12e}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--sizes", type=int, nargs="+", default=[24, 36, 48])
    parser.add_argument("--umax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=400)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()
    scan(args)


if __name__ == "__main__":
    main()
