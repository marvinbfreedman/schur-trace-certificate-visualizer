# Schur Trace Certificate Visualizer

Static GitHub Pages companion archive for the endpoint trace / Schur certificate
visual explainer from the RH Weyl-positive verification work.

Live artifact:

- Visual explainer: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_explainer.html
- Snapshot slideshow: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshot_slideshow.html
- Latest snapshot: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_snapshots/20260628-144016/visual_explainer.html
- Visual data: https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_data.json

Full reproducibility supplement:

- Curated source package: [`supplement/`](supplement/)
- Python quick start: [`supplement/QUICKSTART.md`](supplement/QUICKSTART.md)
- Reproducibility notes: [`supplement/REPRODUCIBILITY.md`](supplement/REPRODUCIBILITY.md)
- Current status: [`supplement/STATUS.md`](supplement/STATUS.md)
- Dashboard report: [`supplement/reports/dashboard.txt`](supplement/reports/dashboard.txt)

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

## Python Quick Start

These commands apply to the full companion source bundle, not only the static
GitHub Pages subset in `docs/`.

Check Python:

```bash
python3 --version
```

Python 3.10 or newer is recommended. Most ledger-reading commands use only the
Python standard library. Numerical scripts may use `numpy`, `scipy`, `mpmath`,
and `matplotlib`.

Install common dependencies when running the full source bundle:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Read the ledger dashboard:

```bash
python3 rh_dashboard.py
python3 rh_dashboard.py --output reports/dashboard.txt --json-out reports/dashboard_summary.json
```

Expected regenerated status from the current archive:

- 293 theorem records are detected.
- 293 theorem records are marked closed.
- 0 theorem-ledger gaps remain open.
- 16 numerical certificates are detected.
- No scan file reports negative inertia/eigenvalue evidence.

Preview which scripts would still run without overwriting existing JSON:

```bash
python3 rh_orchestrator.py --dry-run --skip-existing
```

Run unpaired entry-point scripts with a timeout:

```bash
python3 rh_orchestrator.py --skip-existing --workers 8 --timeout 120
```

Useful filtered runs:

```bash
python3 rh_orchestrator.py --filter scan --dry-run
python3 rh_orchestrator.py --filter certificate --dry-run
python3 rh_orchestrator.py --filter theorem --dry-run
python3 rh_orchestrator.py --filter certificate --workers 8 --timeout 120
python3 rh_orchestrator.py --filter scan --workers 8 --timeout 120
```

The orchestrator sets `PYTHONPATH=source/` and runs each entry-point script in
its own subprocess.

To open the visual explainer locally from a source bundle, use an HTTP server so
the JSON data can be fetched:

```bash
python3 -m http.server 8765
```

Then visit:

```text
http://127.0.0.1:8765/visual/visual_explainer.html
```

## Published Paths

- [`docs/source/visual_explainer.html`](https://marvinbfreedman.github.io/schur-trace-certificate-visualizer/source/visual_explainer.html)
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
