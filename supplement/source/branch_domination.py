#!/usr/bin/env python3
import argparse

from dilation_proof_diagnostics import branch_first_derivative_matrix

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("branch_domination.py requires numpy; run with python") from exc


def eigvals(mat):
    return np.linalg.eigvalsh((mat + mat.T) / 2.0)


def domination_ratios(pos, neg, tol):
    vals, vecs = np.linalg.eigh((neg + neg.T) / 2.0)
    mask = vals < -tol
    if not mask.any():
        return vals, None
    block = vecs[:, mask].T @ ((pos + pos.T) / 2.0) @ vecs[:, mask]
    scale = np.diag(1.0 / np.sqrt(-vals[mask]))
    ratios = np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)
    return vals, ratios


def schur_relative_to_negative_space(total, reference, tol):
    vals, vecs = np.linalg.eigh((reference + reference.T) / 2.0)
    mask = vals < -tol
    if not mask.any():
        return vals, None, None
    e_basis = vecs[:, mask]
    c_basis = vecs[:, ~mask]
    total = (total + total.T) / 2.0
    aa = e_basis.T @ total @ e_basis
    ab = e_basis.T @ total @ c_basis
    bb = c_basis.T @ total @ c_basis
    bb_vals, bb_vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = bb_vals > tol
    if keep.any():
        inv = (bb_vecs[:, keep] / bb_vals[keep]) @ bb_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
    else:
        schur = aa
    range_defect = 0.0
    if (~keep).any():
        null = bb_vecs[:, ~keep]
        range_defect = float(np.linalg.norm(ab @ null))
    return vals, schur, range_defect


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--a", type=float, default=0.245)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=70)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=220)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    q_pos = branch_first_derivative_matrix(
        args.kind, args.a, args.L, args.s_order, args.umax, args.u_order
    )
    q_neg = branch_first_derivative_matrix(
        args.kind, -args.a, args.L, args.s_order, args.umax, args.u_order
    )
    q_sum = 0.5 * (q_pos + q_neg)

    pos_vals = eigvals(q_pos)
    neg_vals, ratios = domination_ratios(q_pos, q_neg, args.tol)
    sum_vals = eigvals(q_sum)
    _, schur, range_defect = schur_relative_to_negative_space(q_sum, q_neg, args.tol)
    print(
        f"branch domination kind={args.kind} a={args.a:g} L={args.L:g} "
        f"s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    print(
        f"  Q_+a: min={pos_vals[0]:.12e} max={pos_vals[-1]:.12e} "
        f"neg={(pos_vals < -args.tol).sum()}"
    )
    print(
        f"  Q_-a: min={neg_vals[0]:.12e} max={neg_vals[-1]:.12e} "
        f"neg={(neg_vals < -args.tol).sum()}"
    )
    print(
        f"  (Q_+a+Q_-a)/2: min={sum_vals[0]:.12e} max={sum_vals[-1]:.12e} "
        f"neg={(sum_vals < -args.tol).sum()}"
    )
    if ratios is not None:
        print(
            f"  Q_+a/(-Q_-a) on Q_-a negative subspace: "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e}"
        )
    if schur is not None:
        schur_vals = eigvals(schur)
        print(
            f"  Schur of combined form over neg(Q_-a): "
            f"min={schur_vals[0]:.12e} max={schur_vals[-1]:.12e} "
            f"range_defect={range_defect:.12e}"
        )


if __name__ == "__main__":
    main()
