#!/usr/bin/env python3
"""Exact derivative diagnostics for b(s)=beta(c e^s)/zeta(c e^s).

The endpoint-zero normalization turns the rows beta,zeta into 1 and

  b(lambda)=beta(lambda)/zeta(lambda),
  beta(lambda)=G(lambda,c)-G(c,c), zeta(lambda)=1-c/lambda.

This script differentiates G(lambda,c) by the J_p recurrence.  It tests whether
b is monotone or Stieltjes-like; that is the natural next analytic hypothesis
for a two-moment quotient proof.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import g_gram, j_asymptotic

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_b_row_derivatives.py requires numpy") from exc


def j_list(a: float, max_half_index: int = 13) -> dict[float, float]:
    """Return J_{k+1/2} for k=0..max_half_index."""
    if a > 80.0:
        return {k + 0.5: j_asymptotic(k + 0.5, a) for k in range(max_half_index + 1)}
    root = math.sqrt(a)
    out = {
        0.5: 1.0 / a
        + 0.5 * math.sqrt(math.pi) * math.exp(a) * math.erfc(root) / (a**1.5)
    }
    for k in range(1, max_half_index + 1):
        p = k + 0.5
        out[p] = 1.0 / a + (p / a) * out[p - 1.0]
    return out


def jp_deriv(js, p: float, order: int) -> float:
    # d^n/da^n J_p = (-1)^n int x^p (x-1)^n exp(-a(x-1)) dx
    #              = sum_{m=0}^n binom(n,m)(-1)^m J_{p+m}.
    total = 0.0
    for m in range(order + 1):
        total += math.comb(order, m) * ((-1.0) ** m) * js[p + m]
    return total


def g_lambda_derivs(lam: float, mu: float, max_order: int = 2):
    """Return d^k/dlambda^k G(lambda,mu), k=0..max_order."""
    a = lam + mu
    js = j_list(a, max_half_index=13)
    vals = []
    for n in range(max_order + 1):
        # G = J3 - a J5 + lam*mu J7.
        # Differentiate by lambda; d/dlambda=d/da on J terms and d lam on lam*mu.
        if n == 0:
            vals.append(js[1.5] - a * js[2.5] + lam * mu * js[3.5])
        elif n == 1:
            vals.append(
                jp_deriv(js, 1.5, 1)
                - js[2.5]
                - a * jp_deriv(js, 2.5, 1)
                + mu * js[3.5]
                + lam * mu * jp_deriv(js, 3.5, 1)
            )
        elif n == 2:
            vals.append(
                jp_deriv(js, 1.5, 2)
                - 2.0 * jp_deriv(js, 2.5, 1)
                - a * jp_deriv(js, 2.5, 2)
                + 2.0 * mu * jp_deriv(js, 3.5, 1)
                + lam * mu * jp_deriv(js, 3.5, 2)
            )
        else:
            raise ValueError("max_order above 2 is not implemented")
    return vals


def b_values(c: float, s: float):
    lam = c * math.exp(s)
    beta = g_gram(lam, c) - g_gram(c, c)
    zeta = 1.0 - c / lam
    g0, g1, g2 = g_lambda_derivs(lam, c, 2)
    beta_lam = g1
    beta_lam2 = g2
    # s derivatives.
    beta_s = lam * beta_lam
    beta_s2 = lam * beta_lam + lam * lam * beta_lam2
    zeta_s = c / lam
    zeta_s2 = -c / lam
    b = beta / zeta
    b_s = (beta_s * zeta - beta * zeta_s) / (zeta * zeta)
    b_s2 = (
        beta_s2 / zeta
        - beta * zeta_s2 / (zeta * zeta)
        - 2.0 * beta_s * zeta_s / (zeta * zeta)
        + 2.0 * beta * zeta_s * zeta_s / (zeta**3)
    )
    d = beta / (lam - c)
    # d = b/lambda because zeta=(lambda-c)/lambda.
    d_s = (b_s - b) / lam
    d_s2 = (b_s2 - 2.0 * b_s + b) / lam
    return b, b_s, b_s2, d, d_s, d_s2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=14.0)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--start", type=float, default=1e-4)
    args = parser.parse_args()

    ss = np.linspace(args.start, args.T, args.samples)
    vals = np.array([b_values(args.c, float(s)) for s in ss])
    b = vals[:, 0]
    b1 = vals[:, 1]
    b2 = vals[:, 2]
    d = vals[:, 3]
    d1 = vals[:, 4]
    d2 = vals[:, 5]
    print(f"endpoint b-row derivatives s in [{args.start:g},{args.T:g}] samples={args.samples}")
    print(f"  b  min={b.min():.12e} max={b.max():.12e}")
    print(f"  b' min={b1.min():.12e} max={b1.max():.12e}")
    print(f"  b'' min={b2.min():.12e} max={b2.max():.12e}")
    print(f"  d=beta/(lambda-c) min={d.min():.12e} max={d.max():.12e}")
    print(f"  d' min={d1.min():.12e} max={d1.max():.12e}")
    print(f"  d'' min={d2.min():.12e} max={d2.max():.12e}")
    print("  selected:")
    for idx in np.linspace(0, args.samples - 1, min(10, args.samples), dtype=int):
        print(
            f"    s={ss[idx]:.8g} b={b[idx]:.12e} "
            f"b'={b1[idx]: .12e} b''={b2[idx]: .12e} "
            f"d={d[idx]:.12e} d'={d1[idx]: .12e}"
        )


if __name__ == "__main__":
    main()
