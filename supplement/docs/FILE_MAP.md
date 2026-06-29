# File Map

## Root

- `rh_dashboard.py`: summarizes JSON theorem/certificate/scan ledgers.
- `rh_orchestrator.py`: runs source scripts in parallel subprocesses.
- `README.md`: repository overview.
- `QUICKSTART.md`: short command guide.
- `REPRODUCIBILITY.md`: reproduction and caveat notes.
- `STATUS.md`: current mathematical status.
- `requirements.txt`: common Python dependencies.
- `CITATION.cff`: citation metadata.

## `source/`

The source directory contains the main Python scripts and JSON ledgers.  Files
are intentionally flat because many scripts were written during an exploratory
research session and import sibling modules directly.

Common naming patterns:

- `*_theorem.py` / `*_theorem.json`: theorem or theorem-wrapper artifacts.
- `*_certificate.py` / `*_certificate.json`: numerical or interval certificate
  artifacts.
- `*_scan.py` / `*_scan.json`: diagnostic scans.
- `*_consequence_theorem.py` / `*_consequence_theorem.json`: narrowed theorem
  interfaces used to reduce dependency-graph risk.

## `reports/`

- `dashboard.txt`: human-readable ledger summary.
- `dashboard_summary.json`: structured dashboard output.

## `visual/`

- `visual_explainer.html`: interactive visual dashboard.
- `visual_data.json`: data loaded by the dashboard.

Run `python3 -m http.server 8765` from the repository root and open
`http://127.0.0.1:8765/visual/visual_explainer.html`.
