#!/usr/bin/env python3
import argparse
import math

from reduced_exact_finite import endpoint_ratios, pieces

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("residual_domination.py requires numpy; run with python, not python3") from exc


def a_value(s, u, pcs):
    es = math.exp(s)
    eu = math.exp(u)
    total = 0.0
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * math.exp(beta * u - c * es * (eu - 1.0))
    return total


def quadrature(end, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return 0.5 * end * (nodes + 1.0), 0.5 * end * weights


def build_matrices(kind, omega, interval_end, s_order, u_end, u_order):
    pcs = pieces(kind)
    s_pts, s_wts = quadrature(interval_end, s_order)
    u_pts, u_wts = quadrature(u_end, u_order)
    root = np.sqrt(s_wts)
    pmat = np.zeros((s_order, s_order), dtype=float)
    rmat = np.zeros_like(pmat)
    hmat = np.zeros_like(pmat)
    for i, s in enumerate(s_pts):
        for j, t in enumerate(s_pts[: i + 1]):
            p_val = 0.0
            r_val = 0.0
            h_val = 0.0
            for u, w in zip(u_pts, u_wts):
                center = float(u + 0.5 * (s + t))
                common = (
                    float(w)
                    * math.cosh(omega * center)
                    * a_value(float(s), float(u), pcs)
                    * a_value(float(t), float(u), pcs)
                )
                p_val += common * float(u + min(s, t))
                r_val += common * float(0.5 * abs(s - t))
                h_val += common
            pmat[i, j] = pmat[j, i] = root[i] * p_val * root[j]
            rmat[i, j] = rmat[j, i] = root[i] * r_val * root[j]
            hmat[i, j] = hmat[j, i] = root[i] * h_val * root[j]
    return pmat, rmat, hmat


def generalized_ratios(pos, residual, tol):
    r_vals, r_vecs = np.linalg.eigh((residual + residual.T) / 2.0)
    neg = r_vals < -tol
    if not neg.any():
        return r_vals, None
    block = r_vecs[:, neg].T @ ((pos + pos.T) / 2.0) @ r_vecs[:, neg]
    scale = np.diag(1.0 / np.sqrt(-r_vals[neg]))
    ratios = np.linalg.eigvalsh((scale @ block @ scale + scale @ block.T @ scale) / 2.0)
    return r_vals, ratios


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--s-order", type=int, nargs="+", default=[40, 60, 80])
    parser.add_argument("--umax", type=float, default=10.0)
    parser.add_argument("--u-order", type=int, default=240)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    print(
        f"residual domination kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} u=[0,{args.umax:g}] u_order={args.u_order}"
    )
    for order in args.s_order:
        pmat, rmat, hmat = build_matrices(
            args.kind, args.omega, args.L, order, args.umax, args.u_order
        )
        full = pmat + rmat
        p_vals = np.linalg.eigvalsh((pmat + pmat.T) / 2.0)
        r_vals, ratios = generalized_ratios(pmat, rmat, args.tol)
        f_vals = np.linalg.eigvalsh((full + full.T) / 2.0)
        h_vals = np.linalg.eigvalsh((hmat + hmat.T) / 2.0)
        print(f"s_order={order}")
        print(
            f"  P_min={p_vals[0]:.12e} R_min={r_vals[0]:.12e} "
            f"Full_min={f_vals[0]:.12e} H_min={h_vals[0]:.12e}"
        )
        print(
            f"  R inertia: neg={(r_vals < -args.tol).sum()} "
            f"zero={(abs(r_vals) <= args.tol).sum()} pos={(r_vals > args.tol).sum()}"
        )
        if ratios is not None:
            print(
                f"  P/(-R) on R-neg: min={ratios[0]:.12e} "
                f"max={ratios[-1]:.12e}"
            )
        print("  R low=", " ".join(f"{v:.4e}" for v in r_vals[:8]))


if __name__ == "__main__":
    main()
