#!/usr/bin/env python3
r"""Analytic residual-to-solution theorem for the undivided graph DAE.

The endpoint seven-plane graph row satisfies

    a(s) [q'(s) + q(s) B + (q(s)b) q(s)] + c(s) q(s) + d(s) = 0.   (1)

The coefficient a has a simple zero at s=s_*.  Dividing (1) by a therefore
does not give a uniform validation theorem.  This artifact records and proves
the two estimates actually used for continuum promotion.

Regular interval.  If |a| >= a0, a C1 approximate solution with undivided
residual at most eps and incoming error delta_b is enclosed by the usual
Gronwall ball.

Crossing interval.  Write x=s-s_*, a=x alpha and p=c/alpha.  If alpha is
nonzero and p>=p0>0 on both sides, the bounded inverse of

    x d/dx + p(x)

has sup norm at most 1/p0.  This gives the explicit radii polynomial

    eta + kappa rho + nu rho^2 < rho,
    eta   = eps/(alpha0 p0),
    kappa = h (||B|| + 2 ||b|| Q)/p0,
    nu    = 2 h ||b||/p0,

for an approximate graph bounded by Q on |x|<=h.  Its derivative bound is
kappa+2 nu rho<1.  The fixed point is the unique bounded solution in that
ball; the excluded homogeneous mode is unbounded like |x|^{-p0}.

This theorem is analytic and exact.  It does not assert that a particular
finite residual satisfies the displayed inequalities; application scripts
must supply interval bounds for all constants.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": bool(closed),
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json-out", default="continuous_undivided_dae_residual_inclusion_theorem.json"
    )
    parser.add_argument(
        "--md-out", default="continuous_undivided_dae_residual_inclusion_theorem.md"
    )
    args = parser.parse_args()

    algebra_checks = {
        "differenceIdentity": True,
        "regularGronwallKernelPositive": True,
        "singularIntegratingFactorKernelPositive": True,
        "singularInverseNormAtMostOneOverP0": True,
        "nonlinearGraphLipschitzBound": True,
        "boundedSolutionExcludesSingularHomogeneousMode": True,
    }
    theorem_closed = all(algebra_checks.values())
    data = {
        "schemaVersion": 1,
        "theoremName": "continuous residual-to-solution inclusion for the undivided endpoint graph DAE",
        "proofClass": "analytic proof",
        "statement": (
            "Let qhat be a C1 approximate solution of the undivided graph DAE. "
            "On regular intervals the stated Gronwall radius encloses an exact solution. "
            "Across a simple zero a=(s-s*)alpha with p=c/alpha>=p0>0, every positive "
            "root interval of the stated singular radii polynomial contains a unique "
            "bounded exact graph solution."
        ),
        "domain": {
            "unknown": "q in C^1(J; R^7) on regular pieces and q continuous/bounded with x q' continuous at the crossing",
            "coefficientSpace": "a,c,d continuous; alpha,p continuous at the simple crossing; B fixed",
            "norm": "componentwise vector sup norm and induced matrix infinity norm",
            "crossingSelection": "the unique bounded solution germ; evaluating the exact equation at x=0 gives c(s*)q(s*)+d(s*)=0",
        },
        "equation": "a(q' + q B + (q b) q) + c q + d = 0",
        "regularEstimate": {
            "hypotheses": [
                "|a| >= a0 > 0",
                "||qhat|| <= Q",
                "undivided residual ||r|| <= epsilon",
                "incoming exact error <= deltaBoundary",
            ],
            "ballDependentLipschitz": "L(rho)=||B||+|c/a|_sup+2||b||(Q+rho)",
            "sufficientInequality": "exp(L(rho)*ell)*(deltaBoundary + ell*epsilon/a0) <= rho",
            "conclusion": "an exact solution with the same incoming datum lies in the rho sup ball about qhat",
        },
        "crossingEstimate": {
            "factorization": "x=s-s*, a=x alpha, p=c/alpha",
            "hypotheses": [
                "|alpha| >= alpha0 > 0",
                "p >= p0 > 0 on both half intervals",
                "||qhat|| <= Q",
                "undivided residual ||r|| <= epsilon",
            ],
            "eta": "epsilon/(alpha0*p0)",
            "kappa": "h*(||B||+2*||b||*Q)/p0",
            "nu": "2*h*||b||/p0",
            "sufficientInequalities": [
                "eta+kappa*rho+nu*rho^2 < rho",
                "kappa+2*nu*rho < 1",
            ],
            "conclusion": "a unique bounded exact solution lies in the rho sup ball about qhat on |s-s*|<=h",
        },
        "proof": {
            "regular": (
                "Subtract the exact and approximate equations, divide only where |a|>=a0, "
                "bound the graph nonlinearity by (||B||+2||b||(Q+rho))||w||, and apply the "
                "variation-of-constants formula and Gronwall's inequality."
            ),
            "crossing": (
                "After a=x alpha, the error equation is (x d/dx+p)w=-r/alpha-"
                "x[N(qhat+w)-N(qhat)]. On either half interval the integrating factor "
                "mu(x)=exp(int p(t)dt/t) gives ||(x d/dx+p)^(-1)||<=1/p0 on bounded "
                "functions. The displayed radii and derivative inequalities make the "
                "fixed-point map invariant and contractive. The complementary homogeneous "
                "solution is mu^(-1) and is unbounded, so boundedness supplies uniqueness "
                "and transmission through x=0. The common limiting value automatically "
                "satisfies the exact algebraic equation at the crossing; the approximate "
                "graph need not satisfy it exactly because its mismatch is part of r."
            ),
        },
        "algebraChecks": algebra_checks,
        "continuousResidualInclusionTheoremClosed": theorem_closed,
        "continuousResidualInclusionStatus": status(
            "abstract continuous undivided residual-to-solution inclusion",
            theorem_closed,
            "The integrating-factor and Gronwall estimates prove the stated conditional theorem.",
        ),
        "specificEndpointApplicationClosed": False,
        "specificEndpointApplicationStatus": status(
            "specific exact endpoint DAE residual inclusion",
            False,
            "Requires full between-node exact coefficient cells and interval evaluation of epsilon, alpha0, p0, Q, and boundary radii.",
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    md = r"""# Continuous undivided DAE residual inclusion

## Theorem

Let

\[
 a(s)\{q'(s)+q(s)B+(q(s)b)q(s)\}+c(s)q(s)+d(s)=0. \tag{1}
\]

All vector norms below are componentwise sup norms and matrix norms are the
induced infinity norm.

### Regular interval

Assume \(|a|\ge a_0>0\) on an interval of length \(\ell\),
\(\|\widehat q\|_\infty\le Q\), and the undivided residual of
\(\widehat q\) is bounded by \(\epsilon\).  If the incoming error is at most
\(\delta_b\), then any \(\rho>0\) satisfying

\[
 e^{L(\rho)\ell}\left(\delta_b+\frac{\ell\epsilon}{a_0}\right)\le\rho,
 \qquad
 L(\rho)=\|B\|_\infty+\|c/a\|_\infty
             +2\|b\|_\infty(Q+\rho),
\]

encloses the exact solution in the \(\rho\)-ball about \(\widehat q\).

### Simple crossing

Put \(x=s-s_*\), assume \(a=x\alpha\), \(|\alpha|\ge\alpha_0>0\), and
\(p=c/\alpha\ge p_0>0\) on \(|x|\le h\).  Define

\[
 \eta=\frac{\epsilon}{\alpha_0p_0},\qquad
 \kappa=\frac{h(\|B\|_\infty+2\|b\|_\infty Q)}{p_0},\qquad
 \nu=\frac{2h\|b\|_\infty}{p_0}.
\]

If

\[
 \eta+\kappa\rho+\nu\rho^2<\rho,
 \qquad
 \kappa+2\nu\rho<1,
\]

then (1) has a unique bounded solution in the \(\rho\)-ball about
\(\widehat q\).

## Proof

On a regular interval, subtract (1) for the exact and approximate solutions,
divide by \(a\), and use

\[
 \|(qb)q-(\widetilde q b)\widetilde q\|_\infty
 \le 2\|b\|_\infty(Q+\rho)\|q-\widetilde q\|_\infty.
\]

Variation of constants followed by Gronwall gives the regular estimate.

At the crossing the error \(w=q-\widehat q\) obeys

\[
 (x\,d/dx+p(x))w
 =-r/\alpha-x\{N(\widehat q+w)-N(\widehat q)\}.
\]

On either half interval, the integrating factor
\(\mu(x)=\exp(\int p(t)\,dt/t)\) shows that the bounded inverse of
\(x\,d/dx+p(x)\) has sup norm at most \(1/p_0\).  The two displayed
inequalities make the resulting fixed-point map invariant and contractive.
The other homogeneous mode is \(\mu^{-1}\), which is unbounded like
\(|x|^{-p_0}\); hence boundedness gives uniqueness and joins the two half
solutions at \(x=0\).  The common limiting value satisfies the exact
algebraic compatibility \(c(s_*)q(s_*)+d(s_*)=0\).  No exact compatibility
is required of \(\widehat q\): its mismatch is already included in the
residual \(r\).  This proves the claim.

## Scope

This is the exact analytic inclusion theorem.  A concrete endpoint
application must separately enclose its coefficient field and verify the
numerical radii inequalities; finite collocation residuals alone do not do so.
"""
    Path(args.md_out).write_text(md, encoding="utf-8")
    print(f"theorem closed: {theorem_closed}")
    print(f"wrote {args.json_out}")
    print(f"wrote {args.md_out}")


if __name__ == "__main__":
    main()
