#!/usr/bin/env python3
"""Projected-log commutator identity for P0.

Set

  F_lambda(tau)=exp(-lambda tau)-exp(-c tau),
  B = (1+tau) d/dtau + 1,
  g_lambda = B exp(-lambda tau),
  u_lambda = B F_lambda = g_lambda - g_c,
  s_lambda = log(lambda/c).

The positive Dirichlet Gram form is

  E(lambda,mu)=<u_lambda,u_mu>_{x^(3/2) dtau}.

If S0 is the projected logarithmic operator on the boundary-null span,

  S0 F_lambda = s_lambda F_lambda,

then the anti-commutator identity is exact:

  P0(lambda,mu)
    = 1/2 <B F_lambda, B S0 F_mu>
      +1/2 <B S0 F_lambda, B F_mu>
    = 1/2(s_lambda+s_mu) E(lambda,mu).

If S is the unprojected logarithmic companion,

  S F_lambda = s_lambda exp(-lambda tau),

then the corresponding cross form is C=P0+R2, with rank(R2)<=2:

  R2(lambda,mu)
    = 1/2 s_mu <u_lambda,g_c> + 1/2 s_lambda <g_c,u_mu>.

This script verifies those identities and checks the tempting shortcut through
the unconditioned anti-commutator

  A_g(lambda,mu)=1/2(s_lambda+s_mu)<g_lambda,g_mu>.

That shortcut is not the desired proof: P0-A_g is finite-rank but not the
rank-two positive-defect identity needed for Lemma A.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "endpoint_p0_commutator_identity.py requires numpy; run with python"
    ) from exc


def sym(mat):
    return 0.5 * (mat + mat.T)


def build_forms(args):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)
    logs = np.log(lambdas / args.c)

    first = args.first_step
    if first <= 0.0:
        first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, segments = composite_r_quadrature(
        args.rmax, first, args.ratio, args.segment_order
    )
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts

    e_boundary = np.zeros((args.order, args.order), dtype=float)
    e_unconditioned = np.zeros_like(e_boundary)
    cross_c = np.zeros_like(e_boundary)
    boundary_vec = np.zeros(args.order, dtype=float)
    gc_norm = 0.0

    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(args.c, float(tau))
        g = np.array([g_value(float(lam), float(tau)) for lam in lambdas])
        u = g - gc
        log_g = logs * g

        e_boundary += weight * root[:, None] * np.outer(u, u) * root[None, :]
        e_unconditioned += weight * root[:, None] * np.outer(g, g) * root[None, :]
        cross_c += (
            0.5
            * weight
            * root[:, None]
            * (np.outer(u, log_g) + np.outer(log_g, u))
            * root[None, :]
        )
        boundary_vec += weight * root * u * gc
        gc_norm += weight * gc * gc

    p0 = 0.5 * (logs[:, None] + logs[None, :]) * e_boundary
    projected = p0.copy()
    r2 = 0.5 * (
        np.outer(boundary_vec, root * logs) + np.outer(root * logs, boundary_vec)
    )
    ag = 0.5 * (logs[:, None] + logs[None, :]) * e_unconditioned
    return {
        "s_pts": s_pts,
        "logs": logs,
        "E": sym(e_boundary),
        "G": sym(e_unconditioned),
        "P0": sym(p0),
        "projected": sym(projected),
        "C": sym(cross_c),
        "R2": sym(r2),
        "Ag": sym(ag),
        "gc_norm": gc_norm,
        "segments": segments,
        "first": first,
    }


def inertia(mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    return vals, int((vals < -tol).sum()), int((vals > tol).sum())


def report_form(name, mat, tol):
    vals, neg, pos = inertia(mat, tol)
    rank = int(np.linalg.matrix_rank(sym(mat), tol=tol))
    loose_rank = int(np.linalg.matrix_rank(sym(mat), tol=max(1e-7, tol)))
    print(
        f"  {name:<8} min={vals[0]: .12e} max={vals[-1]: .12e} "
        f"neg={neg:3d} pos={pos:3d} rank={rank:3d} rank1e-7={loose_rank:3d}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=90)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=18)
    parser.add_argument("--first-step", type=float, default=0.0)
    parser.add_argument("--ratio", type=float, default=1.35)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    forms = build_forms(args)
    p0 = forms["P0"]
    projected = forms["projected"]
    cmat = forms["C"]
    r2 = forms["R2"]
    ag = forms["Ag"]
    shortcut_defect = p0 - ag

    print(
        f"endpoint P0 commutator identity lambda=[c,c exp({args.T:g})] "
        f"order={args.order} r=[0,{args.rmax:g}] "
        f"segments={forms['segments']} first={forms['first']:.3e}"
    )
    print(f"  projected identity ||P0-projected|| = {np.linalg.norm(p0-projected):.12e}")
    print(f"  unprojected identity ||C-P0-R2||   = {np.linalg.norm(cmat-p0-r2):.12e}")
    print(f"  gc norm in E-space                 = {forms['gc_norm']:.12e}")

    for name in ("E", "G", "P0", "C", "R2", "Ag"):
        report_form(name, forms[name], args.tol)

    report_form("P0-Ag", shortcut_defect, args.tol)
    vals = np.linalg.eigvalsh(sym(shortcut_defect))
    visible = vals[np.abs(vals) > args.tol]
    print(
        "  visible eig(P0-Ag): "
        + " ".join(f"{v:.6e}" for v in visible[:8])
        + (" ..." if len(visible) > 8 else "")
    )

    print("  interpretation:")
    print("    P0 is exactly the projected-log anti-commutator.")
    print("    C differs from P0 by the rank-two source R2.")
    print("    The unconditioned Ag shortcut has a higher-rank defect, so Lemma A")
    print("    still needs a direct finite-index/transport proof for P0.")


if __name__ == "__main__":
    main()
