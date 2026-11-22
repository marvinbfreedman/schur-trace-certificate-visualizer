# Augmented-Form Riemann-Kernel Replacement Theorem

This theorem replaces the refuted endpoint-\(B\) positive form.  It does not
rename or silently repair the original claim.

## 1. Coercive replacement form

Let \(T_B\) be the bounded self-adjoint Hilbert--Schmidt operator with kernel
\(K_B\), and put

\[
 \rho=1+\|T_B\|_{\rm HS}.
\]

On \(V_{\rm rep}=H^8(0,L)\), define

\[
 \mathfrak a_{\rm rep}[f]
 =\langle T_Bf,f\rangle
  +\rho\|f\|_2^2
  +\sum_{k=1}^8\|D^kf\|_2^2.                 \tag{1}
\]

Since \(\|T_B\|_{\rm op}\le\|T_B\|_{\rm HS}\),

\[
 \mathfrak a_{\rm rep}[f]
 \ge \|f\|_2^2+\sum_{k=1}^8\|D^kf\|_2^2.
\]

The reverse bound follows from boundedness of \(T_B\).  Thus (1) is a closed
positive form whose norm is equivalent to the standard \(H^8\) norm.  The
certified negative direction of the original form remains present in
\(\langle T_Bf,f\rangle\), but is genuinely dominated by the displayed
regularizer.

## 2. Closed moving trace and high block

For

\[
 \Lambda_a(f)=\sum_{k=0}^8\frac{e_k(a)}{k!}D^kf(a),
\]

pointwise Cauchy--Schwarz and integration over \(I\) give

\[
 \|Rf\|_{L^2(I)}^2
 \le M_e^2\sum_{k=0}^8\|D^kf\|_{L^2(I)}^2
 \le M_e^2\|f\|_{\mathfrak a_{\rm rep}}^2.
\]

Hence \(R:V_{\rm rep}\to L^2(I)\) is bounded and \(\ker R\) is closed.  If
\(F_M\) is any fixed finite-dimensional low/mid space, then

\[
 H_{\rm rep}=\ker R\cap F_M^{\perp_{\mathfrak a_{\rm rep}}}
\]

is a nonvacuous closed Hilbert subspace.  Exterior smooth functions are not
discarded; they are controlled by the augmented energy.

Indeed, every function in
\(C_c^\infty((0.545,L))\) belongs to \(\ker R\).  This is an
infinite-dimensional subspace.  Imposing orthogonality to the fixed
finite-dimensional space \(F_M\) removes at most \(\dim F_M\) directions, so
\(H_{\rm rep}\) remains infinite-dimensional.  Thus the repair does not evade
the counterexample by deleting all exterior-support directions.

## 3. Special source rows

The exact Lagrange source row has the fixed-jet form

\[
 E_uf=
 \left(\sum_{j=0}^7b_j(u)D^jf(s_0),\ p(u)f(s_0)\right).              \tag{2}
\]

One-dimensional Sobolev trace theory makes every jet in (2) bounded on
\(H^8\).  The explicit Volterra/exponential formulas for \(b_j,p\) are
analytic in \(u\); dominated differentiation gives bounded first derivatives
on \(U=[0.08,0.52]\).  Consequently

\[
 \sup_{u\in U}\|E_u\|<\infty,
 \qquad
 \|E_u-E_v\|\le C|u-v|.
\]

Riesz representation therefore gives norm-Lipschitz representers in
\(H_{\rm rep}\).  Moreover, the operator

\[
 E:H_{\rm rep}\longrightarrow L^2(U;\mathbb C^2),
 \qquad (Ef)(u)=E_uf,
\]

has range in the span of the eight functions \((b_j(\cdot),0)\) and the one
function \((0,p(\cdot))\).  It therefore has rank at most nine and is compact
and Hilbert--Schmidt.

## 4. Intrinsic Galerkin/Mosco exhaustion

Choose an orthonormal basis \((\phi_n)\) of \(H_{\rm rep}\), let
\(H_N=\operatorname{span}\{\phi_1,\ldots,\phi_N\}\), and let \(P_N\) be its
orthogonal projection.  The spaces are nested and dense.  They therefore
Mosco-converge to \(H_{\rm rep}\), and \(P_N\to I\) strongly.

Since \(E\) is finite rank,

\[
 \|E(P_N-I)\|\to0,
 \qquad
 \|P_NE^*EP_N-E^*E\|\to0.                         \tag{3}
\]

This is exact source operator-norm consistency; no sampled quadrature is used.

## 5. Exact replacement spectral gap

Put \(S=E^*E\).  If \(S\ne0\), its finite positive spectrum has a smallest
member \(\lambda_* >0\); set \(\delta_{\rm rep}=\lambda_*/2\).  If \(S=0\),
set \(\delta_{\rm rep}=1\).  In both cases \(\delta_{\rm rep}\notin\sigma(S)\).
Equation (3) and the Riesz projection formula give

\[
 \|\mathbf1_{[\delta_{\rm rep},\infty)}(P_NSP_N)
   -\mathbf1_{[\delta_{\rm rep},\infty)}(S)\|\to0.
\]

Thus the replacement continuum \(Q_\delta\) identification closes.

## Exact scope

This is a nonvacuous augmented \(H^8\)-form theorem.  It keeps the original
trace and source formulas and imports the independently closed exact endpoint
right inverse.  It does **not** prove nonnegativity of the original
endpoint-\(B\) quotient form, undo its counterexample, or imply the original
Weyl/KLM or RH conclusion.
