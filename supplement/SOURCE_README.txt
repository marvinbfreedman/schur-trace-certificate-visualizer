RH Weyl-positive companion source bundle
Date: 2026-06-19

Positioning:
This bundle is the companion source archive for an AI-assisted RH verification
system.  The repository is organized as verification infrastructure: Python
theorem generators, certificate builders, scan diagnostics, JSON ledgers,
dashboard summaries, and an orchestrator for reproducible validation runs.

The project used AI as a research accelerator for exploring reductions and
proof routes, but kept the final ledger status grounded in deterministic Python
checks, machine-readable certificates, explicit pass/fail fields, and
reproducible reports.  The workflow is close to evaluation-pipeline thinking:
structured artifacts, provenance, orchestration, dashboards, and auditable
closure criteria over a large artifact graph.

Latest dashboard status:
- The original companion dashboard reported 280/291 theorem records closed, with
  11 open gaps remaining.
- The regenerated reports/dashboard_after_open12.json reports 293/293 theorem
  records closed after additional closure/consequence ledgers were added.
- The archive contains 48,000+ generated source, ledger, certificate, report,
  and snapshot artifacts.

This is a source-only companion bundle for the manuscript
rh_weyl_positive_draft.tex and the submission bundle
submission_bundle_20260618-probe3.

Contents:
- Python verification, scan, certificate, audit, theorem-generation, and
  packaging scripts.
- Theorem and certificate JSON ledgers.
- Publication audit ledgers and dependency graph artifacts.
- TeX manuscript source and existing TeX build metadata/logs.
- Text/HTML regeneration and visual-explainer materials.
- Historical text-only visual snapshot JSON/HTML artifacts.

Excluded intentionally:
- PDFs and generated binary build products.
- Images and screenshots.
- Zip/tar archives.
- Python caches and endpoint cache directories.

This bundle is intended as a companion research/reproducibility archive, not
as the arXiv TeX submission source. The arXiv-facing TeX/PDF pair remains in
submission_bundle_20260618-probe3.
