# RH Weyl/Volterra Companion Supplement

This repository-style supplement accompanies the manuscript
**Finite-core Volterra reductions for a Weyl-positive Riemann phase kernel**.

It is an archival reproducibility package for the reduction, certificate, and
audit ledgers developed in the manuscript.  It is best read as a
reduction/certificate supplement: the files record the Weyl/KLM, Volterra,
trace-frame, Schur, and de Branges bridge structure, together with the scripts
and JSON ledgers used to audit the current proof chain.

The package is meant to expose the proof ledger, certificate records, and
regeneration path directly.  The current regenerated dashboard in this archive
reports all theorem-ledger records closed.

## Quick Start

From the extracted root:

```bash
python3 rh_dashboard.py
python3 rh_dashboard.py --output reports/dashboard.txt --json-out reports/dashboard_summary.json
```

To preview which scripts would still run without overwriting existing JSON:

```bash
python3 rh_orchestrator.py --dry-run --skip-existing
```

To run unpaired entry-point scripts with a timeout:

```bash
python3 rh_orchestrator.py --skip-existing --workers 8 --timeout 120
```

The orchestrator sets `PYTHONPATH=source/` and runs each entry-point script in
its own subprocess.

## Current Dashboard Summary

Generated from this archive:

- JSON files loaded by the current dashboard: 400
- theorem records: 293
- certificate records: 16
- scan/diagnostic records: 66
- other JSON records: 25
- load errors: 0
- theorem records marked closed: 293/293
- open theorem-ledger gaps: 0
- numerical certificates: 16
- scans with negative inertia/eigenvalue evidence: 0

The previous companion dashboard reported 280/291 theorem records closed with
11 open gaps.  The current regenerated dashboard closes the added
closure/consequence ledgers and reports 293/293 closed.

## Repository Layout

```text
.
├── README.md
├── QUICKSTART.md
├── REPRODUCIBILITY.md
├── STATUS.md
├── CITATION.cff
├── requirements.txt
├── rh_dashboard.py
├── rh_orchestrator.py
├── docs/
│   ├── FILE_MAP.md
│   ├── OPEN_GAPS.md
│   └── ZENODO_UPDATE_DESCRIPTION.md
├── reports/
│   ├── dashboard.txt
│   └── dashboard_summary.json
├── source/
│   ├── *.py
│   ├── *.json
│   ├── *.tex
│   └── notes/audit text files
└── visual/
    ├── visual_explainer.html
    └── visual_data.json
```

The heavy generated visual snapshot directory from the earlier raw archive is
intentionally not included here.  The interactive visual explainer and its data
file are retained.

## Main Entry Points

### `rh_dashboard.py`

Reads existing JSON ledgers and summarizes theorem status, certificate status,
scan results, and open gaps.

```bash
python3 rh_dashboard.py --output reports/dashboard.txt --json-out reports/dashboard_summary.json
```

### `rh_orchestrator.py`

Discovers script entry points under `source/`, optionally skips scripts that
already have paired `.json` output, and runs selected scripts in parallel.

```bash
python3 rh_orchestrator.py --dry-run --skip-existing
python3 rh_orchestrator.py --filter certificate --workers 8 --timeout 120
```

## Visual Explainer

Open `visual/visual_explainer.html` from a local HTTP server so that
`visual_data.json` can be fetched:

```bash
python3 -m http.server 8765
```

Then visit:

```text
http://127.0.0.1:8765/visual/visual_explainer.html
```

## Scope Note

This supplement is meant to make the manuscript auditable.  It records what is
proved symbolically, what is interval/certificate backed, what is numerical
evidence, and how the closure ledger can be audited.
