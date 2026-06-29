#!/usr/bin/env python3
r"""Estimate derivative bounds for the moving endpoint trace field.

The trace-resolution lemma will control

    F(a) = Lambda_a(f) = sum_k e_k(a) f^(k)(a)/k!

between sampled points.  Its constants depend on smoothness of the coefficient
field e(a).  This script samples e(a) over the active interval and estimates
sup norms of e, e', and e'' by finite differences.

The output is diagnostic: it identifies the scale of the constants needed in
the trace-resolution proof.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at, dot  # noqa: E402


def f(x):
    return float(x)


def fmt(x, digits=8):
    return mp.nstr(x, digits)


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


def norm2(vec):
    return mp.sqrt(mp.fsum(x * x for x in vec))


def max_abs(vec):
    return max(abs(x) for x in vec) if vec else mp.mpf("0")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s-min", default="0.02")
    parser.add_argument("--s-max", default="0.545")
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--order", type=int, default=50)
    parser.add_argument("--rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="lambda_field_derivative_bounds.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s_min = mp.mpf(args.s_min)
    s_max = mp.mpf(args.s_max)
    if args.samples < 5:
        raise SystemExit("--samples must be at least 5")
    h = (s_max - s_min) / (args.samples - 1)
    centers = [s_min + i * h for i in range(args.samples)]

    vals_rows = []
    vecs = []
    min_gap = mp.inf
    min_cos = mp.mpf("1")
    prev = None
    one_dimensional = True
    for s in centers:
        vals, neg, vec = signed_defect(s, args)
        vals_rows.append(vals)
        vecs.append(vec)
        one_dimensional = one_dimensional and len(neg) == 1
        min_gap = min(min_gap, vals[1] - vals[0])
        if prev is not None:
            min_cos = min(min_cos, abs(dot(prev, vec)))
        prev = vec

    first = []
    second = []
    for i in range(1, args.samples - 1):
        first.append([(vecs[i + 1][k] - vecs[i - 1][k]) / (2 * h) for k in range(args.jet_order)])
        second.append([(vecs[i + 1][k] - 2 * vecs[i][k] + vecs[i - 1][k]) / (h * h) for k in range(args.jet_order)])

    coeff_rows = []
    for k in range(args.jet_order):
        e_vals = [row[k] for row in vecs]
        de_vals = [row[k] for row in first]
        dde_vals = [row[k] for row in second]
        coeff_rows.append(
            {
                "k": k,
                "sup": max_abs(e_vals),
                "infAbs": min(abs(x) for x in e_vals),
                "supD1": max_abs(de_vals),
                "supD2": max_abs(dde_vals),
            }
        )

    e_norm = max(norm2(row) for row in vecs)
    d1_norm = max(norm2(row) for row in first)
    d2_norm = max(norm2(row) for row in second)
    data = {
        "sMin": f(s_min),
        "sMax": f(s_max),
        "samples": args.samples,
        "h": f(h),
        "oneDimensional": one_dimensional,
        "minGap": f(min_gap),
        "minConsecutiveCos": f(min_cos),
        "supNormE": f(e_norm),
        "supNormD1": f(d1_norm),
        "supNormD2": f(d2_norm),
        "coefficients": [
            {
                "k": row["k"],
                "sup": f(row["sup"]),
                "infAbs": f(row["infAbs"]),
                "supD1": f(row["supD1"]),
                "supD2": f(row["supD2"]),
            }
            for row in coeff_rows
        ],
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")

    print(
        f"Lambda derivative scan s=[{s_min},{s_max}] samples={args.samples} "
        f"h={fmt(h)}",
        flush=True,
    )
    print(f"  one_dimensional={one_dimensional}", flush=True)
    print(f"  min_gap={fmt(min_gap, 12)} min_consecutive_cos={fmt(min_cos, 12)}", flush=True)
    print(
        f"  sup ||e||={fmt(e_norm)} sup ||e'||={fmt(d1_norm)} "
        f"sup ||e''||={fmt(d2_norm)}",
        flush=True,
    )
    print("  k   sup|e_k|   inf|e_k|   sup|e_k'|  sup|e_k''|", flush=True)
    for row in coeff_rows:
        print(
            f"  {row['k']:1d} {fmt(row['sup']):>11} {fmt(row['infAbs']):>11} "
            f"{fmt(row['supD1']):>11} {fmt(row['supD2']):>11}",
            flush=True,
        )
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
