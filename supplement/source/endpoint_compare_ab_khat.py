#!/usr/bin/env python3
"""Compare normalized B-branch Khat with the older A-branch endpoint kernel.

The earlier endpoint theorem used A=x d/dtau+3/2.  The current corrected
target uses B=x d/dtau+1 after endpoint normalization:

  Khat_B = Lhat_B + Chat_B.

This script compares Khat_B with the endpoint-normalized A-branch kernel

  Khat_A = int x^(3/2){ log(x)(AFhat)^2 + (AFhat)(AHhat) } dtau.

The exact identity suggested by the numerics is

  Khat_B - Khat_A = 1/2 <Fhat_lambda,Fhat_mu>_{x^(3/2)dtau}.

It follows from A=B+1/2, the endpoint-null integration-by-parts identity
<BF,G>+<F,BG>=-1/2<F,G>, and the commutator B(log(x)F)=log(x)BF+F.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature, feature
from endpoint_c_rank2_split import g_value
from endpoint_closed_form_p0 import beta_closed, p0_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_compare_ab_khat.py requires numpy") from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def v_b(c: float, lam: float, tau: float) -> float:
    return (g_value(lam, tau) - g_value(c, tau)) / (lam - c)


def w_b(c: float, lam: float, tau: float) -> float:
    return math.log(lam / c) * g_value(lam, tau) / (lam - c)


def fhat(c: float, lam: float, tau: float) -> float:
    return (math.exp(-lam * tau) - math.exp(-c * tau)) / (lam - c)


def chat_closed(c: float, lam: float, mu: float) -> float:
    ph = p0_kernel(c, lam, mu) / ((lam - c) * (mu - c))
    d_lam = beta_closed(c, lam) / (lam - c)
    d_mu = beta_closed(c, mu) / (mu - c)
    r_lam = math.log(lam / c) / (lam - c)
    r_mu = math.log(mu / c) / (mu - c)
    return ph + 0.5 * (d_lam * r_mu + r_lam * d_mu)


def build(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    khat_a = np.zeros((args.order, args.order), dtype=float)
    lhat_b = np.zeros_like(khat_a)
    chat_b_quad = np.zeros_like(khat_a)
    fgram = np.zeros_like(khat_a)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        af = np.zeros(args.order, dtype=float)
        ah = np.zeros(args.order, dtype=float)
        vb = np.zeros(args.order, dtype=float)
        wb = np.zeros(args.order, dtype=float)
        for i, lam in enumerate(lambdas):
            afi, ahi, _, _ = feature(args.c, float(lam), float(tau))
            denom = float(lam) - args.c
            af[i] = afi / denom
            ah[i] = ahi / denom
            vb[i] = v_b(args.c, float(lam), float(tau))
            wb[i] = w_b(args.c, float(lam), float(tau))
        fh = np.array([fhat(args.c, float(lam), float(tau)) for lam in lambdas])
        khat_a += (
            weight
            * root[:, None]
            * (
                math.log(x) * np.outer(af, af)
                + 0.5 * (np.outer(af, ah) + np.outer(ah, af))
            )
            * root[None, :]
        )
        lhat_b += (
            weight
            * math.log(x)
            * root[:, None]
            * np.outer(vb, vb)
            * root[None, :]
        )
        chat_b_quad += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(vb, wb) + np.outer(wb, vb))
            * root[None, :]
        )
        fgram += weight * root[:, None] * np.outer(fh, fh) * root[None, :]

    chat_b = np.zeros_like(khat_a)
    for i, lam in enumerate(lambdas):
        for j, mu in enumerate(lambdas[: i + 1]):
            val = root[i] * chat_closed(args.c, float(lam), float(mu)) * root[j]
            chat_b[i, j] = chat_b[j, i] = val
    return (
        sym(khat_a),
        sym(lhat_b + chat_b),
        sym(chat_b_quad - chat_b),
        sym(fgram),
        segments,
        first,
    )


def eig_line(name, mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<18} min={vals[0]: .12e} second={vals[1]: .12e} "
        f"max={vals[-1]: .12e} neg={(vals < -tol).sum():3d}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=8.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--rmax", type=float, default=17.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    ka, kb, chat_defect, fgram, segments, first = build(args)
    print(
        f"endpoint A/B Khat comparison lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] segments={segments} "
        f"first={first:.3e}"
    )
    eig_line("Khat_A", ka, args.tol)
    eig_line("Khat_B", kb, args.tol)
    eig_line("B-A", kb - ka, args.tol)
    eig_line("A-B", ka - kb, args.tol)
    eig_line("1/2 Fgram", 0.5 * fgram, args.tol)
    print(f"  ||Chat_B_quad-closed|| = {np.linalg.norm(chat_defect):.12e}")
    print(f"  ||B-A-1/2Fgram||       = {np.linalg.norm(kb-ka-0.5*fgram):.12e}")


if __name__ == "__main__":
    main()
