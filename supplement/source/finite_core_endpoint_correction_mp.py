#!/usr/bin/env python3
"""Compare full finite-core Volterra jets with the endpoint B-model jet.

This diagnostic computes, in the same formal Taylor-jet basis,

  Delta = K_red(finite core) - K_endpoint,B

and checks whether Delta dominates the negative spectral subspace of the
endpoint B-model.  It is a numerical guide to the analytic correction term
that is absent from the endpoint-only theorem.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_kb_confluent_mp import integrate as endpoint_integrate  # noqa: E402
from tilde3_volterra_confluent_mp import (  # noqa: E402
    integrate as core_integrate,
    pieces_for,
)


def eigvals(mat):
    return mp.eigsy(mp.mpf("0.5") * (mat + mat.T), eigvals_only=True)


def ratios_on_negative(reference, correction):
    vals, vecs = mp.eigsy(mp.mpf("0.5") * (reference + reference.T), eigvals_only=False)
    out = []
    for idx, val in enumerate(vals):
        if val >= 0:
            continue
        v = [vecs[i, idx] for i in range(reference.rows)]
        corr = mp.fsum(
            v[i] * correction[i, j] * v[j]
            for i in range(reference.rows)
            for j in range(reference.cols)
        )
        out.append((idx, val, corr, corr / (-val)))
    return vals, out


def solve_columns(mat, rhs):
    out = mp.matrix(rhs.rows, rhs.cols)
    for col in range(rhs.cols):
        sol = mp.lu_solve(mat, rhs[:, col])
        for row in range(rhs.rows):
            out[row, col] = sol[row]
    return out


def schur_relative_to_reference_negative(reference, total):
    vals, vecs = mp.eigsy(mp.mpf("0.5") * (reference + reference.T), eigvals_only=False)
    neg = [idx for idx, val in enumerate(vals) if val < 0]
    pos = [idx for idx in range(len(vals)) if idx not in neg]
    if not neg or not pos:
        return vals, None, None

    order = neg + pos
    q = mp.matrix(reference.rows, reference.cols)
    for new_col, old_col in enumerate(order):
        for row in range(reference.rows):
            q[row, new_col] = vecs[row, old_col]

    block = q.T * (mp.mpf("0.5") * (total + total.T)) * q
    m = len(neg)
    aa = block[:m, :m]
    ab = block[:m, m:]
    ba = block[m:, :m]
    bb = block[m:, m:]
    bb_vals = eigvals(bb)
    schur = aa - ab * solve_columns(bb, ba)
    schur_vals = eigvals(schur)
    return vals, bb_vals, schur_vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--n", type=int, default=9)
    parser.add_argument("--s0", default="0.5")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--order", type=int, default=70)
    parser.add_argument("--core-umax", default="10")
    parser.add_argument("--endpoint-rmax", default="12")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
    omega = mp.mpf(args.omega)

    endpoint, endpoint_segments = endpoint_integrate(
        "kb", mp.pi, s0, args.n, mp.mpf(args.endpoint_rmax), args.order
    )
    core, core_segments, psi0 = core_integrate(
        s0, omega, pieces_for(args.kind), args.n, mp.mpf(args.core_umax), args.order
    )
    correction = core - endpoint

    endpoint_vals, neg_ratios = ratios_on_negative(endpoint, correction)
    _, core_pos_block_vals, core_schur_vals = schur_relative_to_reference_negative(
        endpoint, core
    )
    core_vals = eigvals(core)
    correction_vals = eigvals(correction)

    print(
        f"finite-core endpoint correction kind={args.kind} omega={omega} "
        f"n={args.n} s0={s0} dps={args.dps} order={args.order}"
    )
    print(
        f"  endpoint_rmax={args.endpoint_rmax} endpoint_segments={endpoint_segments} "
        f"core_umax={args.core_umax} core_segments={core_segments}"
    )
    print(f"  Psi(s0)={mp.nstr(psi0, 25)}")
    print("  endpoint K_B eig low:")
    for val in endpoint_vals[: min(args.n, 5)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  correction Delta eig low:")
    for val in correction_vals[: min(args.n, 5)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  finite-core K_red eig low:")
    for val in core_vals[: min(args.n, 5)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  Delta/(-endpoint) on negative endpoint modes:")
    for idx, val, corr, ratio in neg_ratios:
        print(
            f"    mode={idx} endpoint={mp.nstr(val, 18)} "
            f"Delta={mp.nstr(corr, 18)} ratio={mp.nstr(ratio, 18)}"
        )
    print("  K_red Schur relative to endpoint-negative split:")
    if core_pos_block_vals is None:
        print("    no nontrivial endpoint-negative split")
    else:
        print(f"    positive-block min={mp.nstr(core_pos_block_vals[0], 18)}")
        print(f"    Schur min={mp.nstr(core_schur_vals[0], 18)}")


if __name__ == "__main__":
    main()
