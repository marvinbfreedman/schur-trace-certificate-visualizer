#!/usr/bin/env python3
r"""Probe whether the endpoint Sturm source is a finite s-differential feature.

If the endpoint Sturm Green identity produced a genuine local differential
operator in the spectral variable s, then the source features

    L_log F_s,       L_1 H_s

would lie in the finite span of s-derivatives of the spectral feature pair

    {D_s^j F_s, D_s^j H_s : 0 <= j <= q}

with coefficients depending on s but not on tau.  This is the necessary
algebraic step before extracting a principal coefficient a_q(s).

Analytically this span consists of

    P(tau) exp(-lambda tau) + C exp(-c tau),

where P is a polynomial of degree at most q.  The Sturm sources contain
x^(3/2) and x^(3/2) log(x) factors, so finite membership should fail.  The
script gives a direct least-squares residual for q=10.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_sturm_green_identity import (  # noqa: E402
    f_feature,
    h_feature,
    l_f_feature,
    l_h_feature,
)


def f(x) -> float:
    return float(x)


def fmt(x, digits: int = 8) -> str:
    return mp.nstr(x, digits)


def sample_taus(rmax, count):
    # Log-spaced in r, then converted to tau.  This sees both endpoint and tail.
    if count == 1:
        return [mp.expm1(rmax / 2)]
    return [mp.expm1(rmax * i / (count - 1)) for i in range(count)]


def weighted_lstsq_residual(design, target):
    rows = len(target)
    cols = design.cols
    gram = design.T * design
    rhs = design.T * target
    # Tikhonov only stabilizes the diagnostic; it is far below the displayed
    # residual scale in successful/failing cases.
    trace = mp.fsum(abs(gram[i, i]) for i in range(cols))
    ridge = mp.mpf("1e-60") * max(mp.mpf("1"), trace)
    for i in range(cols):
        gram[i, i] += ridge
    coeff = mp.lu_solve(gram, rhs)
    residual = target - design * coeff
    res_norm = mp.sqrt(mp.fsum(residual[i] ** 2 for i in range(rows)))
    target_norm = mp.sqrt(mp.fsum(target[i] ** 2 for i in range(rows)))
    coeff_norm = mp.sqrt(mp.fsum(coeff[i] ** 2 for i in range(cols)))
    return res_norm / target_norm if target_norm else mp.mpf("0"), coeff_norm


def feature_derivative(which, c, s, tau, order):
    if which == "F":
        return mp.diff(lambda ss: f_feature(c, ss, tau), s, order)
    if which == "H":
        return mp.diff(lambda ss: h_feature(c, ss, tau), s, order)
    raise ValueError(which)


def target_value(which, c, s, tau):
    x = 1 + tau
    if which == "LlogF":
        return l_f_feature(c, s, x, tau, "log")
    if which == "L1H":
        return l_h_feature(c, s, x, tau, "one")
    raise ValueError(which)


def compute_case(args, target_name):
    c = mp.pi
    s = mp.mpf(args.s)
    taus = sample_taus(mp.mpf(args.rmax), args.samples)
    columns = []
    for which in ("F", "H"):
        for order in range(args.max_order + 1):
            columns.append((which, order))

    design = mp.matrix(args.samples, len(columns))
    target = mp.matrix(args.samples, 1)
    for i, tau in enumerate(taus):
        # Mild r-weighting prevents the large-r tail from dominating the
        # least-squares norm completely.
        weight = mp.e ** (-mp.mpf(args.weight_decay) * mp.log1p(tau))
        root = mp.sqrt(weight)
        target[i] = root * target_value(target_name, c, s, tau)
        for j, (which, order) in enumerate(columns):
            design[i, j] = root * feature_derivative(which, c, s, tau, order)

    rel, coeff_norm = weighted_lstsq_residual(design, target)
    return {
        "target": target_name,
        "relativeResidual": f(rel),
        "coefficientNorm": f(coeff_norm),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s", default="0.30")
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--samples", type=int, default=34)
    parser.add_argument("--rmax", default="4")
    parser.add_argument("--weight-decay", default="2")
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--json-out", default="sturm_feature_span_probe.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    rows = [compute_case(args, "L1H"), compute_case(args, "LlogF")]

    print(
        f"Sturm feature-span probe s={args.s} q={args.max_order} "
        f"samples={args.samples}"
    )
    print("  target   relative residual   coeff norm")
    for row in rows:
        print(
            f"  {row['target']:<7} {fmt(mp.mpf(row['relativeResidual']), 10):>18} "
            f"{fmt(mp.mpf(row['coefficientNorm']), 8):>12}"
        )

    data = {
        "s": f(mp.mpf(args.s)),
        "maxOrder": args.max_order,
        "samples": args.samples,
        "rmax": f(mp.mpf(args.rmax)),
        "rows": rows,
        "interpretation": (
            "Large residual means the tau-Sturm source is not generated by a "
            "finite s-differential operator on the feature pair, so no local "
            "principal coefficient a_q(s) can be extracted from that route."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
