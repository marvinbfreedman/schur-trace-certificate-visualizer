#!/usr/bin/env python3
"""Global Galerkin Schur test: endpoint B-model vs finite-core K_red.

This script is the global analogue of finite_core_endpoint_correction_mp.py.
It builds weighted quadrature matrices on s in [0,L] for

  E(s,t) = K_endpoint,B(s,t),
  R(s,t) = K_red(finite core)(s,t),
  Delta = R - E,

and tests R relative to the negative spectral subspace of E.

The endpoint matrix uses the closed mpmath formulas from
endpoint_kb_domination_closed_mp.py.  The finite-core matrix uses the exact
Laguerre formula from reduced_exact_finite.py.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

from endpoint_kb_domination_closed_mp import (  # noqa: E402
    ClosedForms,
    fgram_entry,
    ka_entry,
)
from reduced_exact_finite import full_red, pieces  # noqa: E402


def sym(mat):
    return 0.5 * (mat + mat.T)


def endpoint_kb(cf, c, s, t):
    lam = c * mp.e**s
    mu = c * mp.e**t
    return ka_entry(cf, c, lam, mu) + mp.mpf("0.5") * fgram_entry(cf, c, lam, mu)


def weighted_matrices(kind, omega, length, order, laguerre, dps):
    mp.mp.dps = dps
    nodes, weights = np.polynomial.legendre.leggauss(order)
    pts = 0.5 * length * (nodes + 1.0)
    wts = 0.5 * length * weights
    root = np.sqrt(wts)
    c = mp.pi
    cf = ClosedForms()
    pcs = pieces(kind)
    lag_nodes, lag_weights = np.polynomial.laguerre.laggauss(laguerre)

    endpoint = np.zeros((order, order), dtype=float)
    core = np.zeros_like(endpoint)
    for i, s in enumerate(pts):
        for j, t in enumerate(pts[: i + 1]):
            e = float(endpoint_kb(cf, c, mp.mpf(str(s)), mp.mpf(str(t))))
            r = full_red(float(s), float(t), omega, pcs, lag_nodes, lag_weights)
            endpoint[i, j] = endpoint[j, i] = root[i] * e * root[j]
            core[i, j] = core[j, i] = root[i] * r * root[j]
    return pts, wts, sym(endpoint), sym(core)


def schur_relative_to_negative(reference, total, tol):
    vals, vecs = np.linalg.eigh(sym(reference))
    neg = vals < -tol
    if not neg.any():
        return vals, None, None, None, None
    qn = vecs[:, neg]
    qp = vecs[:, ~neg]
    block = qp.T @ sym(total) @ qp
    block_vals, block_vecs = np.linalg.eigh(sym(block))
    keep = block_vals > tol
    aa = qn.T @ sym(total) @ qn
    ab = qn.T @ sym(total) @ qp
    if keep.any():
        inv = (block_vecs[:, keep] * (1.0 / block_vals[keep])) @ block_vecs[:, keep].T
        schur = aa - ab @ inv @ ab.T
        resid = np.linalg.norm(ab @ (np.eye(len(block_vals)) - block_vecs[:, keep] @ block_vecs[:, keep].T))
    else:
        schur = aa
        resid = np.linalg.norm(ab)
    return vals, block_vals, np.linalg.eigvalsh(sym(schur)), resid, qn


def correction_ratios(reference, correction, tol):
    vals, vecs = np.linalg.eigh(sym(reference))
    out = []
    for idx, val in enumerate(vals):
        if val >= -tol:
            continue
        v = vecs[:, idx]
        corr = float(v.T @ sym(correction) @ v)
        out.append((idx, val, corr, corr / (-val)))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=2.0)
    parser.add_argument("--order", type=int, default=14)
    parser.add_argument("--laguerre", type=int, default=80)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    _, _, endpoint, core = weighted_matrices(
        args.kind, args.omega, args.L, args.order, args.laguerre, args.dps
    )
    delta = core - endpoint
    e_vals, block_vals, schur_vals, resid, _ = schur_relative_to_negative(
        endpoint, core, args.tol
    )
    r_vals = np.linalg.eigvalsh(sym(core))
    d_vals = np.linalg.eigvalsh(sym(delta))
    ratios = correction_ratios(endpoint, delta, args.tol)

    print(
        f"global endpoint/Volterra Schur kind={args.kind} omega={args.omega:g} "
        f"L={args.L:g} order={args.order} laguerre={args.laguerre} dps={args.dps}"
    )
    print(
        f"  endpoint: min={e_vals[0]:.12e} neg={(e_vals < -args.tol).sum()} "
        f"second={e_vals[1]:.12e}"
    )
    print(
        f"  Delta:   min={d_vals[0]:.12e} neg={(d_vals < -args.tol).sum()}"
    )
    print(
        f"  K_red:   min={r_vals[0]:.12e} neg={(r_vals < -args.tol).sum()}"
    )
    print("  Delta/(-endpoint) on endpoint-negative modes:")
    for idx, val, corr, ratio in ratios[:8]:
        print(
            f"    mode={idx} endpoint={val:.12e} "
            f"Delta={corr:.12e} ratio={ratio:.12e}"
        )
    print("  K_red Schur relative to endpoint-negative split:")
    if block_vals is None:
        print("    endpoint has no negative modes at this tolerance")
    else:
        print(f"    positive-block min={block_vals[0]:.12e}")
        print(f"    Schur min={schur_vals[0]:.12e}")
        print(f"    range residual={resid:.12e}")


if __name__ == "__main__":
    main()
