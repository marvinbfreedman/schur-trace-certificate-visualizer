#!/usr/bin/env python3
import argparse

from legendre_certificate import (
    project_to_legendre,
    shifted_legendre_values,
    trial_coefficients,
    weighted_kernel,
)

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("schur_bottleneck.py requires numpy; run with python, not python3") from exc


def complement(q, size):
    full, _ = np.linalg.qr(q, mode="complete")
    return full[:, q.shape[1] :]


def schur_for_tol(kmat, emat, tol):
    rest = complement(emat, kmat.shape[0])
    kk = (kmat + kmat.T) / 2.0
    aa = emat.T @ kk @ emat
    ab = emat.T @ kk @ rest
    bb = rest.T @ kk @ rest
    vals, vecs = np.linalg.eigh((bb + bb.T) / 2.0)
    keep = vals > tol
    if keep.any():
        inv_bb = (vecs[:, keep] * (1.0 / vals[keep])) @ vecs[:, keep].T
        projector = vecs[:, keep] @ vecs[:, keep].T
        sc = aa - ab @ inv_bb @ ab.T
        resid = np.linalg.norm(ab @ (np.eye(len(vals)) - projector))
        min_pos = vals[keep][0]
    else:
        sc = aa
        resid = np.linalg.norm(ab)
        min_pos = 0.0
    sc = (sc + sc.T) / 2.0
    sc_vals, sc_vecs = np.linalg.eigh(sc)
    return {
        "rest": rest,
        "bb_vals": vals,
        "keep_count": int(keep.sum()),
        "min_pos": float(min_pos),
        "range_resid": float(resid),
        "schur": sc,
        "schur_vals": sc_vals,
        "schur_vecs": sc_vecs,
        "aa": aa,
        "ab": ab,
        "inv_bb": inv_bb if keep.any() else np.zeros_like(bb),
    }


def eval_legendre_combo(coeffs, xs, interval_end):
    vals = shifted_legendre_values(xs, len(coeffs), interval_end)
    return vals @ coeffs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--basis", type=int, default=40)
    parser.add_argument("--quad", type=int, default=220)
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument("--trial", choices=("center_exp5", "s_exp1", "one_exp5"), default="center_exp5")
    parser.add_argument("--tols", default="1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15")
    parser.add_argument("--sample-tol", type=float, default=None)
    args = parser.parse_args()

    pts, wts, _, full_w = weighted_kernel(args.kind, args.omega, args.L, args.quad, args.laguerre)
    kmat = project_to_legendre(full_w, pts, wts, args.basis, args.L)
    emat = trial_coefficients(pts, wts, args.basis, args.L, args.trial)
    kvals = np.linalg.eigvalsh((kmat + kmat.T) / 2.0)
    print(
        f"schur bottleneck kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"basis={args.basis} quad={args.quad} laguerre={args.laguerre} trial={args.trial}"
    )
    print(f"K lows: {' '.join(f'{v:.4e}' for v in kvals[:10])}")
    print(f"K highs: {' '.join(f'{v:.4e}' for v in kvals[-6:])}")

    best = None
    for tol_text in args.tols.split(","):
        tol = float(tol_text)
        data = schur_for_tol(kmat, emat, tol)
        sc_vals = data["schur_vals"]
        print(
            f"tol={tol:.0e}: keep={data['keep_count']} min_pos={data['min_pos']:.4e} "
            f"range_resid={data['range_resid']:.4e} "
            f"Schur=({sc_vals[0]:.12e}, {sc_vals[-1]:.12e})"
        )
        if best is None or abs(sc_vals[0]) < abs(best["schur_vals"][0]):
            best = data

    sample_tol = args.sample_tol if args.sample_tol is not None else float(args.tols.split(",")[-1])
    data = schur_for_tol(kmat, emat, sample_tol)
    bb = data["bb_vals"]
    print("BB lows:", " ".join(f"{v:.4e}" for v in bb[:16]))
    print("BB near tol:", " ".join(f"{v:.4e}" for v in bb[(bb > 1e-16) & (bb < 1e-8)][:16]))
    print("AA matrix:")
    for row in data["aa"]:
        print("  " + " ".join(f"{v:.12e}" for v in row))
    print(f"Schur matrix at tol={sample_tol:g}:")
    for row in data["schur"]:
        print("  " + " ".join(f"{v:.12e}" for v in row))

    sc_vec = data["schur_vecs"][:, 0]
    coeff = emat @ sc_vec
    rest_coeff = -data["rest"] @ data["inv_bb"] @ data["ab"].T @ sc_vec
    total_coeff = coeff + rest_coeff
    xs = np.linspace(0.0, args.L, 9)
    samples = eval_legendre_combo(coeff, xs, args.L)
    rest_samples = eval_legendre_combo(rest_coeff, xs, args.L)
    total_samples = eval_legendre_combo(total_coeff, xs, args.L)
    print(
        "small Schur direction in E coefficients="
        + " ".join(f"{v:+.6e}" for v in sc_vec)
    )
    print("E samples:", " ".join(f"({x:.2f},{y:+.6e})" for x, y in zip(xs, samples)))
    print("rest minimizer samples:", " ".join(f"({x:.2f},{y:+.6e})" for x, y in zip(xs, rest_samples)))
    print("total minimizing samples:", " ".join(f"({x:.2f},{y:+.6e})" for x, y in zip(xs, total_samples)))
    print(
        f"L2 coeff norms: E={np.linalg.norm(coeff):.6e} "
        f"rest={np.linalg.norm(rest_coeff):.6e} total={np.linalg.norm(total_coeff):.6e}"
    )


if __name__ == "__main__":
    main()
