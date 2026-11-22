# Current Mathematical Status

Release date: **2026-08-14**

Manuscript: **v2.1.0**

Companion source: **v2.2.0**

## Terminal status

**The Riemann hypothesis is not established by this work.**

The strict root audit is fail-closed. Its Riemann-kernel continuum layer imports
a refuted premise, while three independent layers remain open: the
quotient-to-original Weyl lift, uniform \(\omega\) coverage, and the non-circular
KLM/de Branges bridge.

## Riemann-kernel application

The original endpoint-\(B\) application is **refuted as stated**. A certified
exterior-support direction disproves nonnegativity of the exact Volterra form
on the stated closed-trace core. Two unrelated application obligations also
remain open: compact/Lipschitz continuum source control and exact
Galerkin-to-continuum identification.

The [standalone checker](../docs/source/endpoint-b/verify_endpoint_b.py)
certifies a strictly negative value on 14 exterior nodes. By continuity,
smoothing the witness gives a negative direction with zero trace. On each
nonempty trace fiber that admits this direction, the energy is unbounded
below. [Reproduction instructions](../docs/source/endpoint-b/README.txt)

## Replacement result

The altered-form replacement augments the form with a coercive \(H^8\)
regularizer and closes its ten replacement obligations. This is a theorem about
the reformulated form only. It does not recover the original form, the old
Weyl/KLM chain, or an RH conclusion.

## Zeta transfer and spectral route

The exact identity

\[
W_\zeta(x)=\mathfrak a_{\mathrm{rep}}[Jx]-\mathcal R[Jx]+\mathcal B(x)
\]

has not been closed on exactly matching domains and closures. Negative-direction
reachability and a concrete unconditional zero-location or zero-proportion
conclusion are also unresolved. The limiting prime-side threshold-buffer
theorem remains open; selected finite-height computations are not a uniform
index theorem.

## Controlling artifacts

- `source/strict_rh_prize_proof_audit.json`: `rhConclusionAccepted: false`
- `source/publication_audit_rh_ledger_full.json`: `formalChainClosed: false`
- `source/strict_riemann_kernel_application_audit.json`:
  `riemannKernelApplicationRefutedAsStated: true`
- `source/strict_riemann_kernel_replacement_audit.json`: replacement audit proved
- `source/augmented_to_zeta_transfer_theorem.json`: transfer open
- `source/limiting_prime_side_threshold_buffer_theorem.json`: threshold buffer open

Any historical wrapper, snapshot, or aggregate count that conflicts with these
strict artifacts is non-controlling.
