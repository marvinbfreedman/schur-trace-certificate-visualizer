#!/usr/bin/env python3
"""Closed-form high-precision domination tests for K_B=K_A+1/2 FGram.

This script computes the normalized A-branch endpoint kernel, the positive
divided-difference F-Gram term, and their sum using incomplete-gamma formulas:

  J_p(a)=int_1^inf x^p exp(-a(x-1)) dx
        = exp(a) a^(-p-1) Gamma(p+1,a),

  d_p J_p(a)=int_1^inf x^p log(x) exp(-a(x-1)) dx.

For q_lambda=(3/2-lambda x) exp(-lambda(x-1)),

  Gq(lambda,mu) = int x^(3/2) q_lambda q_mu dx,
  Gqlog(lambda,mu) = int x^(3/2) log(x) q_lambda q_mu dx.

The normalized kernels use F_lambda=(e^{-lambda tau}-e^{-c tau})/(lambda-c).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


class ClosedForms:
    def __init__(self):
        self.j_cache = {}
        self.jlog_cache = {}

    def _key(self, p, a):
        return (str(mp.mpf(p)), mp.nstr(a, mp.mp.dps + 8))

    def j(self, p, a):
        p = mp.mpf(p)
        key = self._key(p, a)
        if key not in self.j_cache:
            self.j_cache[key] = mp.e**a * a ** (-(p + 1)) * mp.gammainc(
                p + 1, a, mp.inf
            )
        return self.j_cache[key]

    def jlog(self, p, a):
        p = mp.mpf(p)
        key = self._key(p, a)
        if key not in self.jlog_cache:
            self.jlog_cache[key] = mp.diff(lambda pp: self.j(pp, a), p)
        return self.jlog_cache[key]

    def gq(self, alpha, beta):
        a = alpha + beta
        return (
            mp.mpf("2.25") * self.j(mp.mpf("1.5"), a)
            - mp.mpf("1.5") * (alpha + beta) * self.j(mp.mpf("2.5"), a)
            + alpha * beta * self.j(mp.mpf("3.5"), a)
        )

    def gqlog(self, alpha, beta):
        a = alpha + beta
        return (
            mp.mpf("2.25") * self.jlog(mp.mpf("1.5"), a)
            - mp.mpf("1.5") * (alpha + beta) * self.jlog(mp.mpf("2.5"), a)
            + alpha * beta * self.jlog(mp.mpf("3.5"), a)
        )

    def fgram_raw(self, alpha, beta):
        a = alpha + beta
        return self.j(mp.mpf("1.5"), a)


def fgram_entry(cf, c, lam, mu):
    return (
        cf.fgram_raw(lam, mu)
        - cf.fgram_raw(lam, c)
        - cf.fgram_raw(c, mu)
        + cf.fgram_raw(c, c)
    ) / ((lam - c) * (mu - c))


def ka_entry(cf, c, lam, mu):
    denom = (lam - c) * (mu - c)
    s_lam = mp.log(lam / c)
    s_mu = mp.log(mu / c)

    lpart = (
        cf.gqlog(lam, mu)
        - cf.gqlog(lam, c)
        - cf.gqlog(c, mu)
        + cf.gqlog(c, c)
    )
    cpart = mp.mpf("0.5") * (
        s_mu * (cf.gq(lam, mu) - cf.gq(c, mu))
        + s_lam * (cf.gq(lam, mu) - cf.gq(lam, c))
    )
    return (lpart + cpart) / denom


def eigvals(mat):
    return mp.eigsy(mat, eigvals_only=True)


def domination_ratio(ka, fg):
    vals, vecs = mp.eigsy(ka, eigvals_only=False)
    neg = [i for i, val in enumerate(vals) if val < 0]
    if not neg:
        return vals, []
    ratios = []
    for idx in neg:
        v = [vecs[j, idx] for j in range(ka.rows)]
        fval = mp.fsum(v[i] * fg[i, j] * v[j] for i in range(ka.rows) for j in range(ka.rows))
        ratios.append((idx, fval / (-2 * vals[idx]), vals[idx], fval))
    return vals, ratios


def relative_threshold(ka, fg):
    """Return eigenvalues of F^(-1/2) K_A F^(-1/2).

    If gamma_min is the first returned eigenvalue, then on this finite span
    K_A + alpha FGram is positive semidefinite exactly when
    alpha >= max(0, -gamma_min), up to the numerical precision of the
    Cholesky solve.
    """

    try:
        chol = mp.cholesky(fg)
    except Exception as exc:  # pragma: no cover - diagnostic path
        return None, f"cholesky failed: {exc}"
    def solve_columns(mat, rhs):
        out = mp.matrix(rhs.rows, rhs.cols)
        for col in range(rhs.cols):
            sol = mp.lu_solve(mat, rhs[:, col])
            for row in range(rhs.rows):
                out[row, col] = sol[row]
        return out

    left = solve_columns(chol, ka)
    rel = solve_columns(chol, left.T).T
    rel = mp.mpf("0.5") * (rel + rel.T)
    return eigvals(rel), None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--s0", type=str, default="0.5")
    parser.add_argument("--h", type=str, default="0.05")
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    c = mp.pi
    s0 = mp.mpf(args.s0)
    h = mp.mpf(args.h)
    lambdas = [c * mp.e ** (s0 + i * h) for i in range(args.n)]
    cf = ClosedForms()
    ka = mp.matrix(args.n)
    fg = mp.matrix(args.n)
    kb = mp.matrix(args.n)
    for i in range(args.n):
        for j in range(i + 1):
            ae = ka_entry(cf, c, lambdas[i], lambdas[j])
            fe = fgram_entry(cf, c, lambdas[i], lambdas[j])
            ka[i, j] = ka[j, i] = ae
            fg[i, j] = fg[j, i] = fe
            kb[i, j] = kb[j, i] = ae + mp.mpf("0.5") * fe

    avals, ratios = domination_ratio(ka, fg)
    rel_vals, rel_error = relative_threshold(ka, fg)
    fvals = eigvals(fg)
    bvals = eigvals(kb)
    print(
        f"endpoint KB closed domination n={args.n} s0={s0} h={h} dps={args.dps}"
    )
    print("  K_A eig low:")
    for val in avals[: min(args.n, 8)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  FGram eig low:")
    for val in fvals[: min(args.n, 5)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  K_B eig low:")
    for val in bvals[: min(args.n, 8)]:
        print(f"    {mp.nstr(val, 25)}")
    print("  FGram/(2*(-K_A)) on negative K_A modes:")
    for idx, ratio, aval, fval in ratios:
        print(
            f"    mode={idx} ratio={mp.nstr(ratio, 25)} "
            f"K_A={mp.nstr(aval, 12)} F={mp.nstr(fval, 12)}"
        )
    print("  generalized K_A relative to FGram:")
    if rel_error is not None:
        print(f"    {rel_error}")
    else:
        gamma = rel_vals[0]
        threshold = max(mp.mpf("0"), -gamma)
        slack = mp.mpf("0.5") - threshold
        print(f"    gamma_min={mp.nstr(gamma, 25)}")
        print(f"    alpha_threshold={mp.nstr(threshold, 25)}")
        print(f"    alpha=1/2 slack={mp.nstr(slack, 25)}")


if __name__ == "__main__":
    main()
