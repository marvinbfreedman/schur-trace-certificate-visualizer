#!/usr/bin/env python3
"""Rational/erfc split for the closed endpoint kernels.

The half-integer recurrence gives

  G(lambda,mu)=<g_lambda,g_mu>
              =G_rat(lambda,mu)+G_erfc(lambda,mu),

where the only non-rational term is R(a)=exp(a) erfc(sqrt(a)),
a=lambda+mu.  This script verifies the split against endpoint_closed_form_p0
and checks whether the finite-rank corrected P0/C kernels become positive
separately on the rational and erfc parts.

The separate positivity check is diagnostic: if it fails, the final Hardy
proof must use cancellation between the rational and erfc pieces.
"""

from __future__ import annotations

import argparse
import math

from endpoint_closed_form_p0 import g_gram

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_closed_form_split.py requires numpy") from exc


def safe_erfc_factor(a: float) -> float:
    """Return exp(a) erfc(sqrt(a)), using the asymptotic branch for large a."""
    if a <= 80.0:
        return math.exp(a) * math.erfc(math.sqrt(a))
    # erfcx(sqrt(a)) ~ 1/(sqrt(pi*a)) sum (-1)^k (2k-1)!!/(2a)^k.
    total = 1.0
    term = 1.0
    for k in range(1, 30):
        term *= -(2 * k - 1) / (2.0 * a)
        total += term
        if abs(term) < 1e-18 * abs(total):
            break
    return total / math.sqrt(math.pi * a)


def g_split(lam: float, mu: float) -> tuple[float, float]:
    a = lam + mu
    lm = lam * mu
    rat = (
        -1.0
        + (lm - 1.5) / a
        + (3.5 * lm - 2.25) / (a * a)
        + 8.75 * lm / (a**3)
        + 13.125 * lm / (a**4)
    )
    erfc = (
        math.sqrt(math.pi)
        * (105.0 * lm - 18.0 * a * a)
        * safe_erfc_factor(a)
        / (16.0 * (a**4.5))
    )
    return rat, erfc


def e_part(c: float, lam: float, mu: float, part: int) -> float:
    def one(x: float, y: float) -> float:
        return g_split(x, y)[part]

    return one(lam, mu) - one(lam, c) - one(c, mu) + one(c, c)


def beta_part(c: float, lam: float, part: int) -> float:
    return g_split(lam, c)[part] - g_split(c, c)[part]


def zeta(c: float, lam: float) -> float:
    return 1.0 - c / lam


def p0_part(c: float, lam: float, mu: float, part: int) -> float:
    return 0.5 * math.log(lam * mu / (c * c)) * e_part(c, lam, mu, part)


def c_part(c: float, lam: float, mu: float, part: int) -> float:
    return p0_part(c, lam, mu, part) + 0.5 * (
        math.log(mu / c) * beta_part(c, lam, part)
        + math.log(lam / c) * beta_part(c, mu, part)
    )


def sym(mat):
    return 0.5 * (mat + mat.T)


def build(args, kernel, beta_kind: str):
    nodes, weights = np.polynomial.legendre.leggauss(args.order)
    s_pts = 0.5 * args.T * (nodes + 1.0)
    s_wts = 0.5 * args.T * weights
    root = np.sqrt(s_wts)
    lambdas = args.c * np.exp(s_pts)

    mats = []
    betas = []
    for part in (0, 1):
        mat = np.zeros((args.order, args.order), dtype=float)
        for i, lam in enumerate(lambdas):
            for j, mu in enumerate(lambdas[: i + 1]):
                val = root[i] * kernel(args.c, float(lam), float(mu), part) * root[j]
                mat[i, j] = mat[j, i] = val
        mats.append(sym(mat))
        if beta_kind == "split":
            betas.append(
                root
                * np.array([beta_part(args.c, float(lam), part) for lam in lambdas])
            )
        else:
            betas.append(np.zeros(args.order, dtype=float))

    full_beta = sum(betas)
    z = root * np.array([zeta(args.c, float(lam)) for lam in lambdas])
    return mats, betas, full_beta, z, lambdas


def eig_report(name, mat, tol):
    vals = np.linalg.eigvalsh(sym(mat))
    print(
        f"  {name:<18} min={vals[0]: .12e} second={vals[1]: .12e} "
        f"neg={(vals < -tol).sum():3d}"
    )
    return vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--T", type=float, default=10.0)
    parser.add_argument("--order", type=int, default=80)
    parser.add_argument("--M", type=float, default=14.0)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    samples = [args.c, args.c * 1.3, args.c * math.exp(2.0), args.c * math.exp(8.0)]
    max_defect = 0.0
    for lam in samples:
        for mu in samples:
            rat, erfc = g_split(lam, mu)
            max_defect = max(max_defect, abs(g_gram(lam, mu) - rat - erfc))

    print(
        f"endpoint rational/erfc split lambda=[c,c exp({args.T:g})] "
        f"order={args.order}"
    )
    print(f"  max pointwise split defect = {max_defect:.12e}")

    for form_name, kernel in (("P0", p0_part), ("C", c_part)):
        mats, betas, full_beta, z, _ = build(args, kernel, "split")
        rat, erfc = mats
        beta_rat, beta_erfc = betas
        full = rat + erfc
        corrected_full = full + args.M * (
            np.outer(full_beta, full_beta) + np.outer(z, z)
        )
        # A natural but non-canonical allocation: beta correction follows the
        # split beta row, while zeta is kept whole with the rational piece.
        corrected_rat = rat + args.M * (
            np.outer(beta_rat, beta_rat) + np.outer(z, z)
        )
        corrected_erfc = erfc + args.M * np.outer(beta_erfc, beta_erfc)
        beta_cross = np.outer(beta_rat, beta_erfc) + np.outer(beta_erfc, beta_rat)
        corrected_erfc_cross = corrected_erfc + args.M * beta_cross

        print(f"  {form_name}:")
        eig_report("rational", rat, args.tol)
        eig_report("erfc", erfc, args.tol)
        eig_report("full", full, args.tol)
        eig_report("corrected full", corrected_full, args.tol)
        eig_report("corrected rat", corrected_rat, args.tol)
        eig_report("corrected erfc", corrected_erfc, args.tol)
        eig_report("erfc + beta cross", corrected_erfc_cross, args.tol)


if __name__ == "__main__":
    main()
