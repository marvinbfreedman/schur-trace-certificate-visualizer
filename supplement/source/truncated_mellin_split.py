#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import psi_mode, quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("truncated_mellin_split.py requires numpy; run with python") from exc


def split_matrices(mode_weights, exponent, interval_end, s_order, u_end, u_order):
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root = np.sqrt(s_wts)
    out_part = np.zeros((s_order, s_order), dtype=float)
    in_part = np.zeros_like(out_part)
    norm = np.zeros_like(out_part)
    for u, w in zip(u_pts, u_wts):
        values = []
        for s in s_pts:
            v = float(s + u)
            total = sum(weight * psi_mode(mode, v) for mode, weight in mode_weights.items())
            values.append(math.exp(exponent * v) * total)
        values = np.array(values, dtype=float)
        gram = np.outer(values, values)
        out_part += float(w * u) * root[:, None] * gram * root[None, :]
        in_part += (
            float(w)
            * root[:, None]
            * (0.5 * (s_pts[:, None] + s_pts[None, :]) * gram)
            * root[None, :]
        )
        norm += float(w) * root[:, None] * gram * root[None, :]
    return out_part, in_part, norm


def eigvals(mat):
    return np.linalg.eigvalsh((mat + mat.T) / 2.0)


def generalized_ratios(pos, neg, tol):
    vals, vecs = np.linalg.eigh((neg + neg.T) / 2.0)
    mask = vals < -tol
    if not mask.any():
        return vals, None
    block = vecs[:, mask].T @ ((pos + pos.T) / 2.0) @ vecs[:, mask]
    scale = np.diag(1.0 / np.sqrt(-vals[mask]))
    ratios = np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)
    return vals, ratios


def schur_relative_to_negative(total, reference, tol):
    vals, vecs = np.linalg.eigh((reference + reference.T) / 2.0)
    mask = vals < -tol
    if not mask.any():
        return None, 0.0
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
    defect = 0.0
    if (~keep).any():
        defect = float(np.linalg.norm(ab @ bb_vecs[:, ~keep]))
    return schur, defect


def parse_weights(text):
    weights = {}
    for part in text.split(","):
        mode_text, weight_text = part.split(":")
        weights[int(mode_text)] = float(weight_text)
    return weights


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", default="1:1")
    parser.add_argument("--a", type=float, default=0.245)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, default=70)
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=220)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    weights = parse_weights(args.weights)
    out_part, in_part, norm = split_matrices(
        weights, args.a, args.L, args.s_order, args.umax, args.u_order
    )
    full = out_part + in_part
    out_vals = eigvals(out_part)
    in_vals, ratios = generalized_ratios(out_part, in_part, args.tol)
    full_vals = eigvals(full)
    norm_vals = eigvals(norm)
    schur, defect = schur_relative_to_negative(full, in_part, args.tol)
    print(
        f"truncated Mellin split weights={weights} a={args.a:g} L={args.L:g} "
        f"s_order={args.s_order} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    print(
        f"  output-log min={out_vals[0]:.12e} max={out_vals[-1]:.12e} "
        f"neg={(out_vals < -args.tol).sum()}"
    )
    print(
        f"  input-log min={in_vals[0]:.12e} max={in_vals[-1]:.12e} "
        f"neg={(in_vals < -args.tol).sum()}"
    )
    print(
        f"  full derivative min={full_vals[0]:.12e} max={full_vals[-1]:.12e} "
        f"neg={(full_vals < -args.tol).sum()}"
    )
    print(f"  norm min={norm_vals[0]:.12e} max={norm_vals[-1]:.12e}")
    if ratios is not None:
        print(
            f"  output-log/(-input-log) on negative input-log: "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e}"
        )
    if schur is not None:
        schur_vals = eigvals(schur)
        print(
            f"  Schur of full over neg(input-log): "
            f"min={schur_vals[0]:.12e} max={schur_vals[-1]:.12e} "
            f"range_defect={defect:.12e}"
        )


if __name__ == "__main__":
    main()
