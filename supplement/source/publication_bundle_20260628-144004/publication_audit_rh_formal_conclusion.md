# Publication Audit Dependency Graph

## Frozen Claim

KLM/Weyl positivity -> augmented closed-cone de Branges kernel positivity -> shifted-Xi endpoint zero exclusion

## Summary

- JSON files scanned: 395
- Theorem/certificate nodes: 314
- Reachable nodes from root: 1
- Edges: 7423
- Formal chain closed in ledger: True
- Independent external proof vetted: False

## Edge Classes

- analytic proof: 3168
- symbolic identity: 2162
- interval/ball certificate: 875
- numerical evidence: 514
- unproven assumption: 704

## Top Blockers


## Danger-Zone Audit

### `klm_debranges_augmented_closed_cone_theorem.json`

- risk: 92
- closed flags: {'statuses.shiftedXiKernelPositivityStatus': True, 'statuses.klmToDeBrangesBridgeStatus': True, 'bridgeClosed': True, 'transformClosed': True}
- break questions:
  - Is the domain exactly the same on both sides?
  - Are closures/density limits valid and explicitly proved?
  - Is positivity strict or only semidefinite, and is strictness ever used?
  - Are boundary terms actually zero in the completed domain?
  - Is any finite certificate promoted to continuum without a compactness theorem?
- top incoming edges:
  - `klm_debranges_bridge_attempt.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 159 (numerical evidence)
  - `klm_debranges_intertwiner_attempt.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 117 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 65 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 65 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 65 (analytic proof)
- top outgoing edges:
  - `klm_debranges_augmented_closed_cone_theorem.json` -> `shifted_xi_debranges_kernel_positivity_theorem.json` risk 77 (symbolic identity)

### `original_weyl_kernel_positivity_assembly_theorem.json`

- risk: 58
- closed flags: {'statuses.branchWeightInputStatus': True, 'statuses.greenContractionInputStatus': True, 'statuses.quadraticFormIdentityInputStatus': True, 'statuses.originalWeylKernelPositivityAssemblyStatus': True, 'originalWeylKernelPositivityAssemblyClosed': True, 'originalWeylKernelPositivityClosed': True, 'originalWeylFormTransportClosed': True, 'greenLiftContractionClosed': True, 'uniformOmegaBranchTransportClosed': True}
- break questions:
  - Is the domain exactly the same on both sides?
  - Are closures/density limits valid and explicitly proved?
  - Is positivity strict or only semidefinite, and is strictness ever used?
  - Are boundary terms actually zero in the completed domain?
  - Is any finite certificate promoted to continuum without a compactness theorem?
- top incoming edges:
  - `publication_audit_external_equivalence_full.json` -> `original_weyl_kernel_positivity_assembly_theorem.json` risk 55 (symbolic identity)
  - `publication_audit_external_equivalence_full.json` -> `original_weyl_kernel_positivity_assembly_theorem.json` risk 55 (symbolic identity)
  - `publication_audit_external_equivalence_full.json` -> `original_weyl_kernel_positivity_assembly_theorem.json` risk 55 (symbolic identity)
  - `publication_audit_external_equivalence_full.json` -> `original_weyl_kernel_positivity_assembly_theorem.json` risk 55 (symbolic identity)
  - `publication_audit_external_equivalence_full.json` -> `original_weyl_kernel_positivity_assembly_theorem.json` risk 55 (symbolic identity)
- top outgoing edges:
  - `original_weyl_kernel_positivity_assembly_theorem.json` -> `green_lift_contraction_consequence_theorem.json` risk 55 (symbolic identity)
  - `original_weyl_kernel_positivity_assembly_theorem.json` -> `original_weyl_quadratic_form_identity_theorem.json` risk 55 (symbolic identity)
  - `original_weyl_kernel_positivity_assembly_theorem.json` -> `original_weyl_branch_weight_theorem.json` risk 50 (interval/ball certificate)

### `xi_augmented_trace_continuum_lift.json`

- risk: 123
- closed flags: {'statuses.muClosedTraceStatus': True, 'statuses.diagonalTailPositiveStatus': True, 'statuses.positivityOnKerRAugStatus': True, 'statuses.finiteAugmentedSchurStatus': True, 'statuses.transportedTraceNormStatus': True, 'statuses.boundedDAugStatus': True, 'statuses.galerkinExhaustionClosureStatus': True, 'continuumAugmentedRepairClosed': True}
- break questions:
  - Is the domain exactly the same on both sides?
  - Are closures/density limits valid and explicitly proved?
  - Is positivity strict or only semidefinite, and is strictness ever used?
  - Are boundary terms actually zero in the completed domain?
  - Is any finite certificate promoted to continuum without a compactness theorem?
- top incoming edges:
  - `klm_debranges_bridge_attempt.json` -> `xi_augmented_trace_continuum_lift.json` risk 109 (analytic proof)
  - `augmented_pullback_limit_theorem.json` -> `xi_augmented_trace_continuum_lift.json` risk 99 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `xi_augmented_trace_continuum_lift.json` risk 77 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `xi_augmented_trace_continuum_lift.json` risk 77 (analytic proof)
  - `publication_audit_external_equivalence_full.json` -> `xi_augmented_trace_continuum_lift.json` risk 77 (analytic proof)
- top outgoing edges:
  - `xi_augmented_trace_continuum_lift.json` -> `xi_augmented_trace_repair_schur.json` risk 109 (analytic proof)
  - `xi_augmented_trace_continuum_lift.json` -> `xi_mellin_convolution_boundary_identity.json` risk 97 (analytic proof)
  - `xi_augmented_trace_continuum_lift.json` -> `xi_mellin_boundary_concomitant.json` risk 97 (analytic proof)
  - `xi_augmented_trace_continuum_lift.json` -> `uniform_omega_weyl_klm_bridge.json` risk 97 (analytic proof)
  - `xi_augmented_trace_continuum_lift.json` -> `continuum_green_lift_closure_theorem.json` risk 82 (interval/ball certificate)

## Downgrade Rules

- Sampling, scans, Galerkin basis evidence, or relative-error checks are numerical evidence unless paired with interval/ball or explicit continuum-closure language.
- Missing JSON imports or open/blocker language become unproven assumptions.
- Closure, density, quotient-domain, endpoint, and boundary-term claims always increase risk even when marked closed.
- Danger-zone nodes and their incident edges receive extra risk weight.
- The audit may report the formal chain closed while still ranking publication blockers.
