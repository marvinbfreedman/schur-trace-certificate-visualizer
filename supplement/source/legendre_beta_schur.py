#!/usr/bin/env python3
import argparse

from legendre_certificate import project_to_legendre, residual_ratio, shifted_legendre_values, weighted_kernel
from schur_bottleneck import schur_for_tol

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("legendre_beta_schur.py requires numpy; run with python, not python3") from exc


def trial_beta(pts, wts, basis_order, interval_end, beta, affine):
    vals = shifted_legendre_values(pts, basis_order, interval_end)
    if affine == "center":
        funcs = [pts - 0.5 * interval_end, np.exp(-beta * pts)]
    elif affine == "one":
        funcs = [np.ones_like(pts), np.exp(-beta * pts)]
    elif affine == "s":
        funcs = [pts.copy(), np.exp(-beta * pts)]
    else:
        raise ValueError(affine)
    coeffs = [vals.T @ (wts * fn) for fn in funcs]
    mat = np.column_stack(coeffs)
    q, r = np.linalg.qr(mat)
    keep = np.abs(np.diag(r)) > 1e-12
    return q[:, keep]


def beta_values(args):
    if args.betas:
        return args.betas
    return np.linspace(args.beta_min, args.beta_max, args.beta_n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.0)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--basis", type=int, default=40)
    parser.add_argument("--quad", type=int, default=220)
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--affine", choices=("center", "one", "s"), default="center")
    parser.add_argument("--beta-min", type=float, default=1.0)
    parser.add_argument("--beta-max", type=float, default=12.0)
    parser.add_argument("--beta-n", type=int, default=45)
    parser.add_argument("--betas", type=float, nargs="*", default=None)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    pts, wts, boundary_w, full_w = weighted_kernel(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    cmat = project_to_legendre(boundary_w, pts, wts, args.basis, args.L)
    kmat = project_to_legendre(full_w, pts, wts, args.basis, args.L)
    tmat = kmat - cmat
    rows = []
    for beta in beta_values(args):
        emat = trial_beta(pts, wts, args.basis, args.L, float(beta), args.affine)
        data = schur_for_tol(kmat, emat, args.tol)
        c_vals, ratios = residual_ratio(cmat, tmat, data["rest"], args.tol)
        data["c_neg"] = int((c_vals < -args.tol).sum())
        data["c_min"] = float(c_vals[0])
        data["ratio_min"] = None if ratios is None else float(ratios[0])
        rows.append((float(data["schur_vals"][0]), float(beta), data))
    rows.sort(reverse=True, key=lambda row: row[0])
    print(
        f"legendre beta Schur kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"basis={args.basis} quad={args.quad} laguerre={args.laguerre} "
        f"tol={args.tol:g} affine={args.affine}"
    )
    for schur_min, beta, data in rows[: args.top]:
        ratio_text = "none" if data["ratio_min"] is None else f"{data['ratio_min']:.3e}"
        print(
            f"beta={beta:.6g} Schur_min={schur_min:.12e} "
            f"Schur_max={data['schur_vals'][-1]:.12e} "
            f"keep={data['keep_count']} min_pos={data['min_pos']:.3e} "
            f"range_resid={data['range_resid']:.3e} "
            f"Cneg={data['c_neg']} Cmin={data['c_min']:.3e} "
            f"tail/res={ratio_text}"
        )


if __name__ == "__main__":
    main()
