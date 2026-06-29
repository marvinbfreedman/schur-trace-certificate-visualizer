# Reproducibility Notes

## What Is Reproducible Here

This supplement contains:

- theorem JSON ledgers;
- numerical certificate JSON ledgers;
- scan/diagnostic JSON ledgers;
- Python scripts that generated or audited those ledgers;
- TeX/manuscript source material retained from the research folder;
- dashboard/orchestrator entry points;
- a lightweight visual explainer.

The package is designed for inspection and regeneration of the certificate
ledger structure, not for a single one-command proof verification.

## Dashboard Reproduction

```bash
python3 rh_dashboard.py --output reports/dashboard.txt --json-out reports/dashboard_summary.json
```

Expected high-level result:

- 400 JSON files load successfully in the current curated supplement;
- 293 theorem records are detected;
- 293 theorem records are marked closed;
- 0 theorem-ledger gaps remain open;
- 16 numerical certificates are detected;
- no scan file reports negative inertia/eigenvalue evidence.

## Script Orchestration

```bash
python3 rh_orchestrator.py --dry-run --skip-existing
python3 rh_orchestrator.py --skip-existing --workers 8 --timeout 120
```

`--skip-existing` leaves already paired `.py`/`.json` artifacts untouched and
runs only scripts that do not already have paired JSON output.

Useful filters:

```bash
python3 rh_orchestrator.py --filter scan --dry-run
python3 rh_orchestrator.py --filter certificate --dry-run
python3 rh_orchestrator.py --filter theorem --dry-run
```

## Numerical Caveat

Some scripts are exploratory or diagnostic.  The theorem ledgers distinguish
between analytic proof, symbolic identity, interval/ball certificate, numerical
evidence, and closure/consequence records.  The dashboard should be treated as
the first source of truth for status.

## Excluded Heavy Files

The earlier raw archive included a large generated `visual_snapshots/` tree.
This repository-style supplement intentionally excludes that directory.  The
interactive visual explainer and `visual_data.json` are included instead.
