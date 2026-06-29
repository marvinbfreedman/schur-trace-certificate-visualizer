#!/usr/bin/env python3
"""Closed-form P0/C kernels from half-integer incomplete-gamma recurrences.

For x=1+tau and

  g_lambda = ((x e^(-lambda tau))') = (1-lambda x)e^(-lambda(x-1)),

define

  G(lambda,mu)=int_1^inf x^(3/2) g_lambda(x) g_mu(x) dx.

Writing a=lambda+mu and

  J_p(a)=int_1^inf x^p e^(-a(x-1)) dx
        =int_0^inf (1+y)^p e^(-a y) dy,

one has

  G(lambda,mu)=J_{3/2}(a)-a J_{5/2}(a)+lambda mu J_{7/2}(a).

The J_p needed here are computed without mpmath.  For moderate a we use the
closed half-integer formula for J_{1/2} and the recurrence

  J_p(a)=1/a + (p/a)J_{p-1}(a).

For large a we use the endpoint asymptotic expansion of J_p.  Since in this
problem a>=2*pi, the moderate branch is well-conditioned.
"""

from __future__ import annotations

import argparse
import math

from endpoint_boundary_null_kernel import composite_r_quadrature
from endpoint_c_rank2_split import g_value

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_closed_form_p0.py requires numpy; run with python") from exc


def binom_real(p: float, k: int) -> float:
    out = 1.0
    for j in range(k):
        out *= (p - j) / (j + 1.0)
    return out


def j_asymptotic(p: float, a: float, terms: int = 28) -> float:
    total = 0.0
    factorial = 1.0
    for k in range(terms):
        if k > 0:
            factorial *= k
        term = binom_real(p, k) * factorial / (a ** (k + 1))
        total += term
        if abs(term) < 1e-18 * max(1.0, abs(total)):
            break
    return total


def j_values(a: float) -> tuple[float, float, float]:
    """Return J_{3/2}, J_{5/2}, J_{7/2}."""
    if a > 80.0:
        return (
            j_asymptotic(1.5, a),
            j_asymptotic(2.5, a),
            j_asymptotic(3.5, a),
        )
    root = math.sqrt(a)
    j_half = 1.0 / a + 0.5 * math.sqrt(math.pi) * math.exp(a) * math.erfc(root) / (
        a ** 1.5
    )
    j_3_2 = 1.0 / a + 1.5 * j_half / a
    j_5_2 = 1.0 / a + 2.5 * j_3_2 / a
    j_7_2 = 1.0 / a + 3.5 * j_5_2 / a
    return j_3_2, j_5_2, j_7_2


def g_gram(lam: float, mu: float) -> float:
    a = lam + mu
    j3, j5, j7 = j_values(a)
    return j3 - a * j5 + lam * mu * j7


def e_kernel(c: float, lam: float, mu: float) -> float:
    return g_gram(lam, mu) - g_gram(lam, c) - g_gram(c, mu) + g_gram(c, c)


def beta_closed(c: float, lam: float) -> float:
    return g_gram(lam, c) - g_gram(c, c)


def zeta_closed(c: float, lam: float) -> float:
    return 1.0 - c / lam


def p0_kernel(c: float, lam: float, mu: float) -> float:
    return 0.5 * math.log(lam * mu / (c * c)) * e_kernel(c, lam, mu)


def c_kernel(c: float, lam: float, mu: float) -> float:
    return p0_kernel(c, lam, mu) + 0.5 * (
        math.log(mu / c) * beta_closed(c, lam)
        + math.log(lam / c) * beta_closed(c, mu)
    )


def sym(mat):
    return 0.5 * (mat + mat.T)


def quadrature_e_beta(c: float, lambdas, rmax: float, segment_order: int):
    first = min(1e-4, 0.05 / float(lambdas[-1]))
    r_pts, r_wts, _ = composite_r_quadrature(rmax, first, 1.35, segment_order)
    tau_pts = np.expm1(r_pts)
    tau_wts = np.exp(r_pts) * r_wts
    n = len(lambdas)
    e_mat = np.zeros((n, n), dtype=float)
    beta = np.zeros(n, dtype=float)
    for tau, tau_wt in zip(tau_pts, tau_wts):
        x = 1.0 + float(tau)
        weight = float(tau_wt) * x**1.5
        gc = g_value(c, float(tau))
        u = np.array(
            [g_value(float(lam), float(tau)) - gc for lam in lambdas],
            dtype=float,
        )
        e_mat += weight * np.outer(u, u)
        beta += weight * u * gc
    return sym(e_mat), beta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--points", type=int, default=12)
    parser.add_argument("--rmax", type=float, default=18.0)
    parser.add_argument("--segment-order", type=int, default=30)
    args = parser.parse_args()

    s_pts = np.linspace(args.T / args.points, args.T, args.points)
    lambdas = args.c * np.exp(s_pts)
    e_quad, beta_quad = quadrature_e_beta(
        args.c, lambdas, args.rmax, args.segment_order
    )
    e_closed = np.array(
        [[e_kernel(args.c, float(lam), float(mu)) for mu in lambdas] for lam in lambdas],
        dtype=float,
    )
    beta = np.array([beta_closed(args.c, float(lam)) for lam in lambdas], dtype=float)
    p0 = np.array(
        [[p0_kernel(args.c, float(lam), float(mu)) for mu in lambdas] for lam in lambdas],
        dtype=float,
    )
    cmat = np.array(
        [[c_kernel(args.c, float(lam), float(mu)) for mu in lambdas] for lam in lambdas],
        dtype=float,
    )
    print(
        f"endpoint closed-form P0/C s=(0,{args.T:g}] points={args.points} "
        f"r=[0,{args.rmax:g}]"
    )
    print(f"  ||E_closed-E_quad||       = {np.linalg.norm(e_closed-e_quad):.12e}")
    print(f"  ||beta_closed-beta_quad|| = {np.linalg.norm(beta-beta_quad):.12e}")
    for name, mat in (("E", e_closed), ("P0", p0), ("C", cmat)):
        vals = np.linalg.eigvalsh(sym(mat))
        print(
            f"  {name:<2} min={vals[0]: .12e} max={vals[-1]: .12e} "
            f"neg={(vals < -1e-10).sum():3d}"
        )
    print("  sample beta/zeta:")
    for s, lam in zip(s_pts[: min(4, len(s_pts))], lambdas[: min(4, len(lambdas))]):
        print(
            f"    s={s:.6g} beta={beta_closed(args.c,float(lam)):.12e} "
            f"zeta={zeta_closed(args.c,float(lam)):.12e}"
        )


if __name__ == "__main__":
    main()
