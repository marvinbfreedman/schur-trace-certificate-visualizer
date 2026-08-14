# Publication Audit Dependency Graph

## Frozen Claim

KLM/Weyl positivity -> augmented closed-cone de Branges kernel positivity -> shifted-Xi endpoint zero exclusion

## Summary

- JSON files scanned: 626
- Theorem/certificate nodes: 539
- Reachable nodes from root: 242
- Edges: 1543
- Formal chain closed in ledger: False
- Formal chain status: refuted-premise
- Strict Riemann-kernel premise refuted as stated: True
- Independent external proof vetted: False

## Edge Classes

- analytic proof: 98
- symbolic identity: 363
- interval/ball certificate: 866
- numerical evidence: 61
- unproven assumption: 155

## Top Blockers

- risk 164: `augmented_to_zeta_transfer_theorem.json` -> `augmented_to_zeta_exact_identity.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, domain/space matching issue, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 144: `asymptotic_zeta_inertia_density_bound.json` -> `asymptotic_zeta_inertia_density_evidence.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, finite-to-continuum promotion issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_002135_0105_q448_d64_r6.json` -> `exact_taylor_6d8db35dd4a4908e9e54f0e8.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_002135_0105_q448_d64_r6.json` -> `exact_taylor_d566a54481053a0b484ebb23.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_00232_widecell_q448_d64_r6.json` -> `exact_taylor_603d3f44e2b1ed3c05d207dc.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0025_widecell_q448_d64_r6.json` -> `exact_taylor_341e093268a964dab7ec7355.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_00288_q448_d64_r6.json` -> `exact_taylor_ac046a34ca029ab876967810.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_00325_narrow_q640_d64_r6.json` -> `exact_taylor_85196c6efbef1f2d53616e0c.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_00345_narrow_q640_d64_r6.json` -> `exact_taylor_33c0a7e44175c19cba57b50e.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0034_q448_d64_r6.json` -> `exact_taylor_c8df1ddbadc6d80660eb7d4f.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0042_00535_q640_d64_r6.json` -> `exact_taylor_6f3792497ae9e547f99d2fd2.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0042_00535_q640_d64_r6.json` -> `exact_taylor_7b8c476a785b8230222223a3.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0042_00535_q640_d64_r6.json` -> `exact_taylor_878ea60291ce55facf5f7a43.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0042_00535_q640_d64_r6.json` -> `exact_taylor_a3223572abd0e60bbd883261.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable
- risk 142: `endpoint_eigenrow_regular_left_0045_q640_d64_r6.json` -> `exact_taylor_4cd4874587a2084e88f879f1.json` (unproven assumption, unknown)
  - reasons: unproven assumption edge, unknown confidence, closure/density/continuum passage issue, boundary term issue
  - classification: referenced JSON is missing or unreadable

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
  - `strict_rh_prize_proof_audit.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 73 (symbolic identity)
  - `strict_klm_debranges_bridge_audit.json` -> `klm_debranges_augmented_closed_cone_theorem.json` risk 63 (symbolic identity)
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
- The audit may report the formal chain closed only when the strict Riemann-kernel application audit is affirmatively closed and not refuted.
