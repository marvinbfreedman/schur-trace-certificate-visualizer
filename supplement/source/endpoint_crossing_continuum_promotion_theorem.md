# Exact endpoint crossing promotion

Let (J(s)) be the exact confluent endpoint matrix and let (e(s)) be its
normalized negative eigenrow.  On

\[
 I_*=[0.15125,,0.191808480511711849\ldots],
\]

the degree-40 Arb certificate supplies eight overlapping complex disks.  An
analytic Newton--Krawczyk argument closes on every disk.  If (h) is the real
cell half-width, (r) the inner complex radius, and \(\rho_j\) the analytic
eigenpair radius on disk (j), Cauchy's estimate gives

\[
 \|e^{(q)}-\widehat e_j^{(q)}\|_{\infty,I_j}
 \leq \frac{q!\rho_j}{(r-h)^q},\qquad 0\leq q\leq8.
\]

The closed cells overlap and cover (I_*), so these are genuine
between-node enclosures, not node samples.

The derivative cells determine interval coefficients (a,c,d) for the exact
undivided graph equation

\[
 a\{q'+qB+(qb)q\}+cq+d=0.
\]

Writing (a=(s-s_*)\alpha) and (p=c/\alpha), the concrete enclosure gives

\[
 |\alpha|\geq 2.6433348820\times10^{-5},\qquad
 p\geq7.2741714715.
\]

With Taylor scaling radius (0.05), the exact residual is at most
(3.000250997\times10^{-6}).  The singular fixed-point constants are

\[
 \eta=0.01560349433,\quad
 \kappa=0.6563411237,\quad
 \nu=1.047510988.
\]

The selected radius (ho=0.1088732726) satisfies

\[
 \eta+\kappa\rho+\nu\rho^2<\rho,
 \qquad \kappa+2\nu\rho=0.8844330222<1.
\]

Therefore the exact undivided DAE has a unique bounded solution in that tube
through the simple crossing.  This theorem closes the continuum promotion on
(I_*).  It does not, by itself, propagate the plane across the two remaining
regular intervals or certify the final restricted active endpoint minor.
