#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("dy_hardy_split.py requires numpy; run with python") from exc


def b_value(p, c, x):
    return math.exp(p * x - c * math.exp(x))


def b_derivative(p, c, x):
    return (p - c * math.exp(x)) * b_value(p, c, x)


def sym(mat):
    return (mat + mat.T) / 2.0


def eigvals(mat):
    return np.linalg.eigvalsh(sym(mat))


def generalized_ratios(pos, neg, tol):
    vals, vecs = np.linalg.eigh(sym(neg))
    mask = vals < -tol
    if not mask.any():
        return vals, None
    block = vecs[:, mask].T @ sym(pos) @ vecs[:, mask]
    scale = np.diag(1.0 / np.sqrt(-vals[mask]))
    ratios = np.linalg.eigvalsh(sym(scale @ block @ scale))
    return vals, ratios


def schur_relative_to_negative(total, reference, tol):
    vals, vecs = np.linalg.eigh(sym(reference))
    mask = vals < -tol
    if not mask.any():
        return None, 0.0, 0
    e_basis = vecs[:, mask]
    c_basis = vecs[:, ~mask]
    total = sym(total)
    aa = e_basis.T @ total @ e_basis
    ab = e_basis.T @ total @ c_basis
    bb = c_basis.T @ total @ c_basis
    bb_vals, bb_vecs = np.linalg.eigh(sym(bb))
    keep = bb_vals > tol
    if keep.any():
        inv = (bb_vecs[:, keep] / bb_vals[keep]) @ bb_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
    else:
        schur = aa
    defect = 0.0
    if (~keep).any():
        defect = float(np.linalg.norm(ab @ bb_vecs[:, ~keep]))
    return schur, defect, int(keep.sum())


def build_matrices(p, c, input_end, input_order, output_end, output_order):
    s_pts, s_wts = quadrature(input_end, input_order)
    r_pts, r_wts = quadrature(output_end, output_order)
    root = np.sqrt(s_wts)
    s_sum = 0.5 * (s_pts[:, None] + s_pts[None, :])
    p_dy = np.zeros((input_order, input_order), dtype=float)
    r_dy_deta = np.zeros_like(p_dy)
    p_plus_r = np.zeros_like(p_dy)
    y_l_y = np.zeros_like(p_dy)
    y_eta = np.zeros_like(p_dy)
    y_norm = np.zeros_like(p_dy)

    # Endpoint products at r=0.
    b0 = np.array([b_value(p, c, float(s)) for s in s_pts])
    endpoint_y_eta = root[:, None] * (0.5 * (s_pts[:, None] + s_pts[None, :])) * np.outer(
        b0, b0
    ) * root[None, :]

    for r, w in zip(r_pts, r_wts):
        b = np.array([b_value(p, c, float(r + s)) for s in s_pts])
        db = np.array([b_derivative(p, c, float(r + s)) for s in s_pts])
        bb = root[:, None] * np.outer(b, b) * root[None, :]
        ddb = root[:, None] * np.outer(db, db) * root[None, :]

        p_dy += float(w * r) * ddb
        r_dy_deta += float(w) * s_sum * ddb
        p_plus_r += float(w) * (float(r) + s_sum) * ddb
        y_l_y += float(w * r) * bb
        y_eta += float(w) * s_sum * bb
        y_norm += float(w) * bb

    a_const = 2.0 * p - 3.0
    lower = (
        a_const * a_const * (y_l_y + y_eta)
        + 2.0 * a_const * y_norm
        + 2.0 * a_const * endpoint_y_eta
    )
    full = lower + 4.0 * p_plus_r
    return {
        "P_DyLDy": p_dy,
        "R_DyDeta": r_dy_deta,
        "P_plus_R": p_plus_r,
        "A2_y": a_const * a_const * (y_l_y + y_eta),
        "2A_norm": 2.0 * a_const * y_norm,
        "2A_endpoint": 2.0 * a_const * endpoint_y_eta,
        "lower": lower,
        "full_expanded": full,
    }


def report_matrix(name, mat, tol):
    vals = eigvals(mat)
    print(
        f"  {name:<14} min={vals[0]: .12e} max={vals[-1]: .12e} "
        f"neg={(vals < -tol).sum():3d}"
    )
    return vals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=1.495)
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--L", type=float, default=12.0)
    parser.add_argument("--s-order", type=int, default=80)
    parser.add_argument("--umax", type=float, default=12.0)
    parser.add_argument("--u-order", type=int, default=260)
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--show-modes", type=int, default=0)
    args = parser.parse_args()

    s_pts, s_wts = quadrature(args.L, args.s_order)
    mats = build_matrices(
        args.p, args.c, args.L, args.s_order, args.umax, args.u_order
    )
    print(
        f"Dy Hardy split p={args.p:g} c={args.c:g} input=[0,{args.L:g}] "
        f"s_order={args.s_order} output=[0,{args.umax:g}] u_order={args.u_order}"
    )
    for name in (
        "P_DyLDy",
        "R_DyDeta",
        "P_plus_R",
        "A2_y",
        "2A_norm",
        "2A_endpoint",
        "lower",
        "full_expanded",
    ):
        report_matrix(name, mats[name], args.tol)

    _, ratios = generalized_ratios(mats["P_DyLDy"], mats["R_DyDeta"], args.tol)
    if ratios is not None:
        print(
            "  P_DyLDy/(-R_DyDeta) on neg(R): "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e}"
        )

    _, ratios = generalized_ratios(4.0 * mats["P_plus_R"], mats["lower"], args.tol)
    if ratios is not None:
        print(
            "  4(P+R)/(-lower) on neg(lower): "
            f"min={ratios[0]:.12e} max={ratios[-1]:.12e}"
        )

    schur, defect, rank = schur_relative_to_negative(
        mats["full_expanded"], mats["lower"], args.tol
    )
    if schur is not None:
        vals = eigvals(schur)
        print(
            "  Schur full over neg(lower): "
            f"min={vals[0]:.12e} max={vals[-1]:.12e} "
            f"bb_rank={rank} range_defect={defect:.12e}"
        )

    if args.show_modes:
        vals, vecs = np.linalg.eigh(sym(mats["lower"]))
        root = np.sqrt(s_wts)
        print("  negative lower modes:")
        for j in range(min(args.show_modes, len(vals))):
            if vals[j] >= -args.tol:
                break
            coeff = vecs[:, j] / root
            norm = math.sqrt(float(np.sum(s_wts * coeff * coeff)))
            coeff = coeff / norm
            moments = [
                float(np.sum(s_wts * (s_pts ** k) * coeff)) for k in range(4)
            ]
            endpoint = float(coeff[0])
            endpoint_slope = float((coeff[1] - coeff[0]) / (s_pts[1] - s_pts[0]))
            print(
                f"    mode {j}: eig={vals[j]:.12e} "
                f"F(0~)={endpoint:.12e} F'(0~)={endpoint_slope:.12e} "
                f"moments0..3={' '.join(f'{m:.6e}' for m in moments)}"
            )


if __name__ == "__main__":
    main()
