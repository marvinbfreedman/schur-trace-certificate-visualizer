# Schur Trace Certificate Visualizer

Static GitHub Pages companion archive for the endpoint trace / Schur certificate
visual explainer from the RH Weyl-positive verification work.

Live artifact:

- Visual explainer: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_explainer.html
- Snapshot slideshow: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshot_slideshow.html
- Latest snapshot: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshots/20260628-144016/visual_explainer.html
- Visual data: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_data.json

## RH Weyl-Positive Companion Source Bundle

Date: 2026-06-19

This repository publishes the static visualizer and machine-readable source
artifacts for an AI-assisted RH verification system. The underlying companion
bundle is organized as verification infrastructure: Python theorem generators,
certificate builders, scan diagnostics, JSON ledgers, dashboard summaries, and
an orchestrator for reproducible validation runs.

The project used AI as a research accelerator for exploring reductions and proof
routes, but kept the final ledger status grounded in deterministic Python checks,
machine-readable certificates, explicit pass/fail fields, and reproducible
reports. The workflow is close to evaluation-pipeline thinking: structured
artifacts, provenance, orchestration, dashboards, and auditable closure criteria
over a large artifact graph.

## Latest Dashboard Status

- The original companion dashboard reported 280/291 theorem records closed, with
  11 open gaps remaining.
- The regenerated `reports/dashboard_after_open12.json` reports 293/293 theorem
  records closed after additional closure/consequence ledgers were added.
- The archive contains 48,000+ generated source, ledger, certificate, report,
  and snapshot artifacts in the source bundle.

This is a source-only companion publication for the manuscript
`rh_weyl_positive_draft.tex` and the frozen submission bundle
`submission_bundle_20260618-probe3`.

## Published Paths

- `docs/source/visual_explainer.html`
- `docs/source/visual_snapshot_slideshow.html`
- `docs/source/visual_data.json`
- `docs/source/rounded_ball_krawczyk_enclosure.svg`
- `docs/source/visual_snapshots/`

## Contents

- Static HTML visual explainer and snapshot slideshow.
- Theorem and certificate JSON ledgers used by the visualizer.
- Publication audit ledgers and dependency graph artifacts.
- Text/HTML regeneration and visual-explainer materials.
- Historical text-only visual snapshot JSON/HTML artifacts.

## Excluded Intentionally

- PDFs and generated binary build products.
- Images and screenshots outside the static visualizer assets.
- Zip/tar archives.
- Python caches and endpoint cache directories.

This repo intentionally serves raw static output from `docs/` with
`docs/.nojekyll`. It is a companion research/reproducibility archive, not the
arXiv TeX submission source. The arXiv-facing TeX/PDF pair remains in
`submission_bundle_20260618-probe3`.
