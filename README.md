# Schur Trace Certificate Visualizer

Static GitHub Pages visualizer and curated audit mirror for the manuscript
**Finite-core Volterra reductions for a Weyl-positive Riemann phase kernel**.

## Current releases

- Manuscript v2.1.0: [doi:10.5281/zenodo.21936429](https://doi.org/10.5281/zenodo.21936429)
- Companion source archive v2.2.0: [doi:10.5281/zenodo.21936544](https://doi.org/10.5281/zenodo.21936544)
- Manuscript concept DOI: [doi:10.5281/zenodo.20746137](https://doi.org/10.5281/zenodo.20746137)
- Source-archive concept DOI: [doi:10.5281/zenodo.20757731](https://doi.org/10.5281/zenodo.20757731)

The Zenodo source archive is the canonical frozen release. This repository
keeps the static visualizer and a small, review-critical subset of the current
audit artifacts; it intentionally does not duplicate the full archive.

## Current mathematical status

The v2.1.0/v2.2.0 release is explicitly fail-closed:

- the strict RH root audit reports **not established**;
- the original Riemann-kernel application is **refuted as stated** because its
  closed-trace nonnegativity premise is disproved;
- the coercive \(H^8\)-augmented form closes only the altered-form replacement
  audit and does not reinstate the original Weyl/KLM or RH claim;
- the exact augmented-to-zeta transfer and the limiting spectral
  threshold-buffer route remain **open**.

Accordingly, this repository and the associated releases do **not** claim a
proof of the Riemann hypothesis. Historical snapshots and aggregate dashboard
counts are non-controlling; the strict root and application audits control.

Reviewer starting points:

- [Standalone 14-node endpoint-B checker](docs/source/endpoint-b/verify_endpoint_b.py) and [reproduction instructions](docs/source/endpoint-b/README.txt)
- [`supplement/STATUS.md`](supplement/STATUS.md)
- [`strict_rh_prize_proof_audit.json`](supplement/source/strict_rh_prize_proof_audit.json)
- [`strict_riemann_kernel_application_audit.json`](supplement/source/strict_riemann_kernel_application_audit.json)
- [`augmented_to_zeta_transfer_theorem.json`](supplement/source/augmented_to_zeta_transfer_theorem.json)
- [`limiting_prime_side_threshold_buffer_theorem.json`](supplement/source/limiting_prime_side_threshold_buffer_theorem.json)

## Visualizer

- [Visual explainer](https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_explainer.html)
- [Snapshot slideshow](https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshot_slideshow.html)
- [June 28 historical snapshot](https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshots/20260628-144016/visual_explainer.html)
- [Visual data](https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_data.json)

The visualizer presents finite certificate and diagnostic material. It is not a
terminal proof audit and must be read together with the current fail-closed
ledgers above.

## Repository layout

- `docs/`: GitHub Pages output and current controlling JSON mirrors.
- `supplement/`: curated human-readable status and audit mirrors.
- Zenodo v2.2.0: complete source-only reproducibility archive, including Python,
  theorem ledgers, deterministic outputs, manuscript TeX, and Lean sources.

Release date: 2026-08-14.
