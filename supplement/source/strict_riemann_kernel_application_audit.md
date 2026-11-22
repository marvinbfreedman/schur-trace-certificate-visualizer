# Strict Riemann-Kernel Application Audit

**Status:** `refuted-as-stated`

The four abstract mechanisms are available, but not all application hypotheses have been proved for the exact Riemann-kernel form.

| Application obligation | Status | Missing proof |
|---|---|---|
| closability of the exact endpoint-B Volterra form | proved | Closed by exact_volterra_form_closability_theorem.md. |
| nonnegativity of the exact Volterra form on the closed-trace core | disproved | The present obligation cannot be proved: reformulate the form, trace, or domain to exclude the certified exterior negative direction. |
| exact endpoint fundamental flow and active singular-value bound | proved | Closed by endpoint_global_dae_right_inverse_continuum_certificate.json. |
| compact and Lipschitz continuum source operator in H_A | open | Prove the special source-row estimate ||E_u f||^2<=eta(u)a[f] on ker R with locally bounded eta, and prove the corresponding representers are norm-Lipschitz in u. |
| Galerkin-to-continuum identification of Q_delta | open | Prove Mosco convergence of the actual closed-trace high blocks and operator-norm convergence of the compact source compressions; then apply Riesz-projection continuity across a certified spectral gap. |

## Wrapper consistency

The strict result follows `high_block_compact_exhaustion_proof.json`. Any broader high-block wrapper is accepted only when it also imports the primary compact-source and completed-form source estimates.
