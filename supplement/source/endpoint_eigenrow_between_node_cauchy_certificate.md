# Between-node Bernstein/Cauchy eigenrow certificate

Let (J(s)) be the exact (9\times9) confluent endpoint matrix and let
((e(s),\lambda(s))) denote its normalized negative eigenpair.  Put

\[
 F(s,e,\lambda)=
 \begin{pmatrix}
 J(s)e-\lambda e\\[2mm]
 (e^{\mathsf T}e-1)/2
 \end{pmatrix}.
\]

At a real center (s_0), the certificate encloses the Taylor coefficients of
(J) through degree (N) by Arb quadrature.  The quadrature remainder is
bounded on explicit Bernstein ellipses and the (r)-integral tail is bounded
analytically.  Near (s=0) it uses the cancellation-preserving identities

\[
 v=e^{-c\tau}\left[-(1-cx)\tau\,
   \operatorname{exprel}(-\delta\tau)-x e^{-\delta\tau}\right],
 \qquad
 w=\frac{s}{c(e^s-1)}g,
 \]

where \(c=\pi\), \(x=e^r\), \(\tau=x-1\), and
\(\delta=c(e^s-1)\).  Both right-hand sides are analytic at \(s=0\); this
avoids the false estimate obtained by bounding a vanishing numerator and
(c(e^s-1)) separately.

Let ((\widehat e,\widehat\lambda)) be the degree-(N) Taylor candidate and
let (H=D_{(e,\lambda)}F).  A polynomial matrix (R(s)) is obtained by
recursively inverting the Taylor series of (H(s,\widehat e,
\widehat\lambda)).  On (|s-s_0|\le r), set

\[
 Y=\|R F(s,\widehat e,\widehat\lambda)\|_\infty,
 \qquad
 Z=\|I-RH(s,\widehat e,\widehat\lambda)\|_\infty.
\]

The finite polynomial parts are evaluated with Arb.  The omitted Taylor tail
of (J) is bounded from an exact majorant (M_R) on a larger disk of radius
(R>r):

\[
 \left\|J-\sum_{k=0}^{N}J_k(s-s_0)^k\right\|_{\max}
 \le
 4M_R\frac{(r/R)^{N+1}}{1-r/R}.
\]

For the bordered eigenpair system, changing the unknown by at most
(\rho) changes the Jacobian by at most (9\rho) in the induced infinity
norm.  Hence the analytic Newton map is invariant and contractive whenever

\[
 Y+Z\rho+9\|R\|_\infty\rho^2<\rho,
 \qquad
 Z+9\|R\|_\infty\rho<1.
\]

This proves the existence and uniqueness of the normalized analytic
eigenpair in the displayed tube.  If the certified real cell has half-width
(h<r), Cauchy's estimate gives, for (0\le q\le8),

\[
 \left\|e^{(q)}-\widehat e^{(q)}\right\|_\infty
 \le \frac{q!\rho}{(r-h)^q}.
\]

The JSON output records these componentwise derivative intervals.  A global
between-node conclusion is emitted only if the union of closed real cells
covers the entire requested interval.  Local disk closure alone is not
reported as global continuation.
