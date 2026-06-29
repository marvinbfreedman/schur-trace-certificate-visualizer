#!/usr/bin/env python3
"""Track the endpoint B-model Taylor-jet defect family.

For each center s0, this script builds the endpoint B-model confluent matrix

  J_n(s0) = [d_s^i d_t^j K_endpoint,B(s0,s0)/(i! j!)]_{i,j=0}^{n-1},

extracts its negative eigendirection, aligns signs across s0, and reports the
Taylor-normalized functional

  Lambda_s0(f) = sum_k e_k(s0) f^(k)(s0)/k!.

The goal is to identify whether the endpoint defect is a one-dimensional
smooth family of distributional/confluent functionals.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_kb_confluent_mp import integrate  # noqa: E402


def dot(a, b):
    return mp.fsum(x * y for x, y in zip(a, b))


def norm(a):
    return mp.sqrt(dot(a, a))


def normalize(vec):
    nrm = norm(vec)
    return [x / nrm for x in vec]


def defect_at(s0, n, rmax, order, tol):
    mat, _ = integrate("kb", mp.pi, s0, n, rmax, order)
    vals, vecs = mp.eigsy(mat, eigvals_only=False)
    neg = [idx for idx, val in enumerate(vals) if val < -tol]
    idx = 0
    vec = normalize([vecs[row, idx] for row in range(n)])
    return vals, neg, vec


def format_vec(vec, digits):
    return "[" + ", ".join(mp.nstr(v, digits) for v in vec) + "]"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=9)
    parser.add_argument("--s-min", default="0.35")
    parser.add_argument("--s-max", default="0.65")
    parser.add_argument("--samples", type=int, default=13)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--order", type=int, default=70)
    parser.add_argument("--rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--digits", type=int, default=12)
    parser.add_argument("--show-coeffs", action="store_true")
    parser.add_argument("--root-iters", type=int, default=18)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s_min = mp.mpf(args.s_min)
    s_max = mp.mpf(args.s_max)
    rmax = mp.mpf(args.rmax)
    tol = mp.mpf(args.tol)
    if args.samples == 1:
        centers = [s_min]
    else:
        step = (s_max - s_min) / (args.samples - 1)
        centers = [s_min + i * step for i in range(args.samples)]

    rows = []
    prev = None
    max_angle = mp.mpf("0")
    min_abs_cos = mp.mpf("1")
    max_coeff_step = mp.mpf("0")
    one_dimensional = True
    min_gap = mp.inf

    for s0 in centers:
        vals, neg, vec = defect_at(s0, args.n, rmax, args.order, tol)
        if prev is not None and dot(prev, vec) < 0:
            vec = [-x for x in vec]
        if prev is not None:
            cos = max(mp.mpf("-1"), min(mp.mpf("1"), dot(prev, vec)))
            abs_cos = abs(cos)
            min_abs_cos = min(min_abs_cos, abs_cos)
            max_angle = max(max_angle, mp.acos(abs_cos))
            max_coeff_step = max(
                max_coeff_step, max(abs(vec[i] - prev[i]) for i in range(args.n))
            )
        prev = vec
        one_dimensional = one_dimensional and len(neg) == 1
        gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
        min_gap = min(min_gap, gap)
        rows.append((s0, vals, neg, vec))

    roots = []
    for (left_s, left_vals, _, _), (right_s, right_vals, _, _) in zip(
        rows[:-1], rows[1:]
    ):
        left = left_vals[0]
        right = right_vals[0]
        if left == 0:
            roots.append(left_s)
            continue
        if left * right > 0:
            continue
        lo_s, hi_s = left_s, right_s
        lo_v, hi_v = left, right
        for _ in range(args.root_iters):
            mid_s = (lo_s + hi_s) / 2
            mid_vals, _, _ = defect_at(mid_s, args.n, rmax, args.order, tol)
            mid_v = mid_vals[0]
            if lo_v * mid_v <= 0:
                hi_s, hi_v = mid_s, mid_v
            else:
                lo_s, lo_v = mid_s, mid_v
        roots.append((lo_s + hi_s) / 2)

    print(
        f"endpoint defect family n={args.n} s=[{s_min},{s_max}] "
        f"samples={args.samples} dps={args.dps} order={args.order} rmax={rmax}"
    )
    print(f"  tol={tol}")
    print(f"  one_dimensional={one_dimensional}")
    print(f"  min_gap(lambda_1-lambda_0)={mp.nstr(min_gap, 20)}")
    print(f"  min_abs_consecutive_cos={mp.nstr(min_abs_cos, 20)}")
    print(f"  max_consecutive_angle={mp.nstr(max_angle, 20)}")
    print(f"  max_coeff_step={mp.nstr(max_coeff_step, 20)}")
    if roots:
        print("  lambda0_zero_crossings=" + ", ".join(mp.nstr(root, 20) for root in roots))
    else:
        print("  lambda0_zero_crossings=none on scanned interval")
    print("  rows:")
    for idx, (s0, vals, neg, vec) in enumerate(rows):
        gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
        print(
            f"    s0={mp.nstr(s0, 8)} neg={len(neg)} "
            f"lambda0={mp.nstr(vals[0], 18)} "
            f"lambda1={mp.nstr(vals[1], 18)} "
            f"gap={mp.nstr(gap, 12)}"
        )
        if args.show_coeffs or idx in (0, len(rows) // 2, len(rows) - 1):
            print(f"      e={format_vec(vec, args.digits)}")

    print("  Lambda_s0 convention:")
    print("    Lambda_s0(f) = sum_k e_k(s0) f^(k)(s0)/k! with ||e||_2=1")


if __name__ == "__main__":
    main()
