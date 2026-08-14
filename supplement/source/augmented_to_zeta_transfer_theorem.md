# Exact Augmented-to-Zeta Transfer Gate

## Status

This theorem is **open**.  The coercive augmented form is valid, but no zeta
conclusion may be inferred from its positivity until every item below is
proved on exactly matching domains and closures.

## Required identity

Let \(X_\zeta\) be a precisely specified zeta/Weil test space and let

\[
 J:X_\zeta\longrightarrow H^8(0,L)
\]

be the exact pullback.  The required identity is

\[
 W_\zeta(x)
 =\mathfrak a_{\rm rep}[Jx]-\mathcal R[Jx]+\mathcal B(x),
\]

where

\[
 \mathcal R[f]
 =(1+\|K_B\|_{\rm HS})\|f\|_2^2
  +\sum_{k=1}^8\|D^k f\|_2^2.
\]

The term \(\mathcal B\) must contain every endpoint, primitive-prefix,
completion, and theta-tail contribution.  The identity must hold first on an
explicit common core and then pass to the stated closures in a proved common
topology.  Finite pullback identities, entrywise kernel limits, and density of
a transported trace image do not substitute for this theorem.

## Reachability fork

The certified endpoint-\(B\) negative direction must be tested against
\(\overline{J(X_\zeta)}^{H^8}\).

- If it belongs to that closure, the positivity route is refuted on the zeta
  image.
- If it is separated from that closure, the separation and the resulting
  natural restricted-domain theorem must be exact.
- Without either result the transfer remains open.

The present counterexample certificate proves existence of a smooth
exterior-support negative witness by mollification, but it does not yet export
one explicit smooth witness with an \(H^8\) enclosure.  The existing pullback
artifacts do not characterize the full \(H^8\) range closure.  Reachability is
therefore currently unresolved.

## Quantitative target

Because the original form is signed, the preferred target is a finite-height
negative-index bound followed by a rigorous asymptotic passage to a specific
unconditional zero-location or zero-proportion theorem.  No internal form
closure is accepted as the terminal conclusion.

## Fail-closed rule

`augmented_to_zeta_transfer_theorem.json` is `proved` only if the exact
pullback/domain identity, regularizer subtraction, boundary/tail accounting,
reachability fork, quantitative inertia theorem, and concrete zeta conclusion
all close.  Missing or malformed inputs leave the theorem `open`.
