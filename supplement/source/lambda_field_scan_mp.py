#!/usr/bin/env python3
"""Scan the endpoint-defect jet hyperplane field.

The continuum constraint

  Lambda_a(f) = sum_{k=0}^8 e_k(a) f^(k)(a)/k! = 0

is a rank-one condition on the 9-jet of f at every active a.  A tempting
shortcut is to solve it as an eighth-order ODE using the coefficient e_8(a),
but this script checks whether that coefficient is a valid global pivot.

The sign convention is e_0(a) < 0.  This is stable in all active scans and
removes the arbitrary eigenvector sign.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at, dot  # noqa: E402


def signed_defect(s, args):
    vals, neg, vec = defect_at(
        s,
        args.jet_order,
        mp.mpf(args.rmax),
        args.order,
        mp.mpf(args.tol),
    )
    if vec[0] > 0:
        vec = [-x for x in vec]
    return vals, neg, vec


def parse_coeffs(text):
    if not text:
        return []
    return [int(piece) for piece in text.replace(",", " ").split()]


def bisect_zero(k, left_s, right_s, left_value, right_value, args):
    f_left = left_value
    f_right = right_value
    if f_left == 0:
        return left_s
    if f_right == 0:
        return right_s
    if f_left * f_right > 0:
        return None
    lo, hi = left_s, right_s
    flo, fhi = f_left, f_right
    for _ in range(args.root_iters):
        mid = (lo + hi) / 2
        _, _, v_mid = signed_defect(mid, args)
        fmid = v_mid[k]
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return (lo + hi) / 2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s-min", default="0.02")
    parser.add_argument("--s-max", default="0.545")
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--order", type=int, default=55)
    parser.add_argument("--rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--digits", type=int, default=12)
    parser.add_argument("--zero-coeffs", default="8")
    parser.add_argument("--root-iters", type=int, default=18)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s_min = mp.mpf(args.s_min)
    s_max = mp.mpf(args.s_max)
    if args.samples == 1:
        centers = [s_min]
    else:
        step = (s_max - s_min) / (args.samples - 1)
        centers = [s_min + i * step for i in range(args.samples)]

    rows = []
    prev = None
    min_abs = [mp.inf for _ in range(args.jet_order)]
    max_abs = [mp.mpf("0") for _ in range(args.jet_order)]
    min_abs_cos = mp.mpf("1")
    one_dimensional = True
    min_gap = mp.inf

    for s in centers:
        vals, neg, vec = signed_defect(s, args)
        if prev is not None:
            min_abs_cos = min(min_abs_cos, abs(dot(prev, vec)))
        prev = vec
        gap = vals[1] - vals[0]
        min_gap = min(min_gap, gap)
        one_dimensional = one_dimensional and len(neg) == 1
        for k, coeff in enumerate(vec):
            min_abs[k] = min(min_abs[k], abs(coeff))
            max_abs[k] = max(max_abs[k], abs(coeff))
        pivot = max(range(args.jet_order), key=lambda k: abs(vec[k]))
        rows.append((s, vals, neg, vec, pivot))

    zero_coeffs = parse_coeffs(args.zero_coeffs)
    zero_roots = {k: [] for k in zero_coeffs}
    for k in zero_coeffs:
        for left, right in zip(rows[:-1], rows[1:]):
            left_value = left[3][k]
            right_value = right[3][k]
            if left_value * right_value > 0:
                continue
            root = bisect_zero(k, left[0], right[0], left_value, right_value, args)
            if root is not None:
                zero_roots[k].append(root)

    print(
        f"Lambda field scan s=[{s_min},{s_max}] samples={args.samples} "
        f"dps={args.dps} jet_order={args.jet_order} order={args.order}"
    )
    print("  sign convention: e_0(a) < 0")
    print(f"  one_dimensional={one_dimensional}")
    print(f"  min_gap(lambda_1-lambda_0)={mp.nstr(min_gap, 20)}")
    print(f"  min_abs_consecutive_cos={mp.nstr(min_abs_cos, 20)}")
    print("  coefficient absolute ranges:")
    for k in range(args.jet_order):
        print(
            f"    k={k}: min_abs={mp.nstr(min_abs[k], 16)} "
            f"max_abs={mp.nstr(max_abs[k], 16)}"
        )
    print("  zero crossings:")
    for k in zero_coeffs:
        roots = zero_roots[k]
        if roots:
            print(
                f"    e_{k}=0 at "
                + ", ".join(mp.nstr(root, 20) for root in roots)
            )
        else:
            print(f"    e_{k}=0: none detected")
    print("  rows:")
    for s, vals, neg, vec, pivot in rows:
        coeffs = ", ".join(mp.nstr(x, args.digits) for x in vec)
        print(
            f"    s={mp.nstr(s, 10)} neg={len(neg)} "
            f"lambda0={mp.nstr(vals[0], 16)} pivot=e_{pivot} "
            f"e=[{coeffs}]"
        )


if __name__ == "__main__":
    main()
