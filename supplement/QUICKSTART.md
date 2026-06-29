# Quick Start

These commands assume you are in the extracted supplement root.

## 1. Check Python

```bash
python3 --version
```

Python 3.10 or newer is recommended.

## 2. Install Common Dependencies

Most ledger-reading commands only need the Python standard library.  Numerical
scripts may use `numpy`, `scipy`, `mpmath`, and `matplotlib`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## 3. Read the Ledger Dashboard

```bash
python3 rh_dashboard.py
```

To write files:

```bash
python3 rh_dashboard.py --output reports/dashboard.txt --json-out reports/dashboard_summary.json
```

## 4. Preview Remaining Script Runs

```bash
python3 rh_orchestrator.py --dry-run --skip-existing
```

This reports entry-point scripts without paired JSON output.

## 5. Run a Filtered Batch

Examples:

```bash
python3 rh_orchestrator.py --filter certificate --workers 8 --timeout 120
python3 rh_orchestrator.py --filter scan --workers 8 --timeout 120
```

## 6. Open the Visual Explainer

Use a local server rather than opening the HTML file directly:

```bash
python3 -m http.server 8765
```

Then open:

```text
http://127.0.0.1:8765/visual/visual_explainer.html
```
