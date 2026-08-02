# Continuous undivided DAE residual inclusion

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
