#!/usr/bin/env python3
"""
rh_dashboard.py — Summarize all RH verification certificate ledgers.

Reads every .json file in the source directory, classifies it by schema type,
and prints a structured status report.

Usage:
    python rh_dashboard.py [--source-dir PATH] [--output FILE] [--json]
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path


# ── Classification ────────────────────────────────────────────────────────────

def classify(name: str, data: dict) -> str:
    if "_load_error" in data:
        return "load_error"
    keys = set(data.keys())
    if "theoremName" in keys:
        return "theorem"
    if "positiveModes" in keys and "rank" in keys:
        return "certificate"
    if "rank" in keys and "nullity" in keys:
        return "certificate"
    if "rows" in keys or "gamma2Top" in keys or "kMin" in keys or "inertia" in keys:
        return "scan"
    return "other"


# ── Theorem analysis ──────────────────────────────────────────────────────────

def analyze_theorems(items):
    closed = []
    open_ = []
    proof_classes = Counter()
    status_counts = Counter()

    for name, d in items:
        gap = d.get("remainingAnalyticGap")
        pc = d.get("proofClass") or "unknown"
        proof_classes[pc] += 1

        statuses = d.get("statuses") or {}
        for key, val in statuses.items():
            if isinstance(val, dict):
                s = val.get("status", "unknown")
                status_counts[s] += 1

        if gap:
            open_.append({"file": name, "gap": gap, "proofClass": pc})
        else:
            closed.append(name)

    return {
        "total": len(items),
        "closed": len(closed),
        "open": len(open_),
        "proof_classes": dict(proof_classes.most_common()),
        "status_counts": dict(status_counts.most_common()),
        "open_items": open_,
    }


# ── Certificate analysis ───────────────────────────────────────────────────────

def analyze_certificates(items):
    rows_out = []
    for name, d in items:
        rows_out.append({
            "file": name,
            "basis": d.get("basis"),
            "rank": d.get("rank"),
            "nullity": d.get("nullity"),
            "positiveModes": d.get("positiveModes"),
            "constraints": d.get("constraints"),
            "evalRangeRelativeDefect": d.get("evalRangeRelativeDefect"),
            "maxFactorRelativeDefect": d.get("maxFactorRelativeDefect"),
        })
    return rows_out


# ── Scan/diagnostic analysis ──────────────────────────────────────────────────

def analyze_scans(items):
    rows_out = []
    neg_any = 0

    for name, d in items:
        rows = d.get("rows") or []
        kmin_vals = []
        neg_counts = []
        pos_counts = []
        gamma_vals = []

        for row in rows:
            if not isinstance(row, dict):
                continue
            kmin = row.get("kMin")
            if kmin is not None:
                kmin_vals.append(kmin)
            kermin = row.get("kerMin")
            if kermin is not None:
                kmin_vals.append(kermin)
            inertia = row.get("inertia")
            if isinstance(inertia, dict):
                neg_counts.append(inertia.get("neg", 0))
                pos_counts.append(inertia.get("pos", 0))
            if "positiveModes" in row:
                pos_counts.append(row["positiveModes"])
            g = row.get("gamma2Top") or row.get("gamma2")
            if g is not None:
                gamma_vals.append(g)

        entry = {
            "file": name,
            "num_rows": len(rows),
            "gap": d.get("gap"),
            "status": d.get("status"),
        }
        if kmin_vals:
            entry["kMin_min"] = min(kmin_vals)
            entry["kMin_max"] = max(kmin_vals)
        if neg_counts:
            entry["max_neg_inertia"] = max(neg_counts)
            entry["all_psd"] = max(neg_counts) == 0
            if max(neg_counts) > 0:
                neg_any += 1
        if gamma_vals:
            entry["gamma2_max"] = max(abs(g) for g in gamma_vals)

        rows_out.append(entry)

    return rows_out, neg_any


# ── Formatting ────────────────────────────────────────────────────────────────

def fmt_sci(val):
    if val is None:
        return "—"
    try:
        return f"{val:.2e}"
    except Exception:
        return str(val)


def fmt_bool(val):
    if val is None:
        return "—"
    return "YES" if val else "NO"


def print_report(source_dir: Path, groups: dict, out):
    theorems = analyze_theorems(groups["theorem"])
    certs = analyze_certificates(groups["certificate"])
    scans, neg_any = analyze_scans(groups["scan"])
    other = groups["other"]
    errors = groups["load_error"]

    sep = "=" * 72

    def p(*args):
        print(*args, file=out)

    p(sep)
    p("RH WEYL-POSITIVE VERIFICATION — CERTIFICATE DASHBOARD")
    p(f"Source: {source_dir}")
    p(sep)
    p()

    # ── Overall ──────────────────────────────────────────────────────────────
    total_json = sum(len(g) for g in groups.values())
    p("OVERALL LEDGER")
    p("-" * 40)
    p(f"  Total JSON files:    {total_json}")
    p(f"  Theorem records:     {len(groups['theorem'])}")
    p(f"  Certificate records: {len(groups['certificate'])}")
    p(f"  Scan/diagnostic:     {len(groups['scan'])}")
    p(f"  Other:               {len(groups['other'])}")
    p(f"  Load errors:         {len(errors)}")
    p()

    # ── Theorems ─────────────────────────────────────────────────────────────
    p("THEOREMS")
    p("-" * 40)
    p(f"  Total:  {theorems['total']}")
    p(f"  Closed: {theorems['closed']}  ({100*theorems['closed']/max(theorems['total'],1):.1f}%)")
    p(f"  Open:   {theorems['open']}")
    p()
    p("  Proof classes:")
    for cls, cnt in sorted(theorems["proof_classes"].items(), key=lambda x: -x[1]):
        p(f"    {cnt:5d}  {cls}")
    p()

    if theorems["open_items"]:
        p("  Open gaps:")
        for item in theorems["open_items"]:
            gap_short = item["gap"][:100].replace("\n", " ")
            p(f"    [{item['proofClass']}] {item['file']}")
            p(f"      Gap: {gap_short}")
        p()

    # ── Certificates ──────────────────────────────────────────────────────────
    p("NUMERICAL CERTIFICATES")
    p("-" * 40)
    if not certs:
        p("  (none)")
    else:
        header = f"  {'File':<52} {'rank':>4} {'null':>4} {'pos':>4} {'evalDefect':>12} {'factorDefect':>13}"
        p(header)
        p("  " + "-" * 90)
        for c in sorted(certs, key=lambda x: x["file"]):
            name = c["file"].replace("_certificate.json", "").replace(".json", "")
            p(
                f"  {name:<52} {str(c['rank'] or '—'):>4} "
                f"{str(c['nullity'] or '—'):>4} "
                f"{str(c['positiveModes'] or '—'):>4} "
                f"{fmt_sci(c['evalRangeRelativeDefect']):>12} "
                f"{fmt_sci(c['maxFactorRelativeDefect']):>13}"
            )
    p()

    # ── Scans ─────────────────────────────────────────────────────────────────
    p("SCAN / DIAGNOSTIC RESULTS")
    p("-" * 40)
    p(f"  Total scan files: {len(scans)}")
    p(f"  Files with negative inertia: {neg_any}")
    p(f"  Files with explicit 'gap' field set: "
      f"{sum(1 for s in scans if s.get('gap'))}")
    p()

    # Show scans that have kMin info (Weyl-positivity direct checks)
    kmin_scans = [s for s in scans if "kMin_min" in s]
    if kmin_scans:
        p(f"  Kernel-minimum scans ({len(kmin_scans)} files):")
        p(f"  {'File':<55} {'kMin_min':>14} {'kMin_max':>14}")
        p("  " + "-" * 85)
        for s in sorted(kmin_scans, key=lambda x: x["file"]):
            name = s["file"].replace(".json", "")
            p(
                f"  {name:<55} "
                f"{fmt_sci(s['kMin_min']):>14} "
                f"{fmt_sci(s['kMin_max']):>14}"
            )
        p()

    # Scans with inertia data
    inertia_scans = [s for s in scans if "max_neg_inertia" in s]
    if inertia_scans:
        bad = [s for s in inertia_scans if not s.get("all_psd")]
        p(f"  Inertia-checked scans: {len(inertia_scans)}  (PSD failures: {len(bad)})")
        if bad:
            p("  NON-PSD SCANS:")
            for s in bad:
                p(f"    {s['file']}  max_neg={s['max_neg_inertia']}")
        p()

    # ── Other ─────────────────────────────────────────────────────────────────
    if other:
        p("OTHER JSON FILES")
        p("-" * 40)
        for name, d in other:
            p(f"  {name}  keys: {sorted(d.keys())[:6]}")
        p()

    # ── Errors ────────────────────────────────────────────────────────────────
    if errors:
        p("LOAD ERRORS")
        p("-" * 40)
        for name, d in errors:
            p(f"  {name}: {d.get('_load_error', '?')}")
        p()

    p(sep)
    p(f"Summary: {theorems['closed']}/{theorems['total']} theorems closed, "
      f"{len(certs)} certificates, "
      f"{neg_any} scans with negative eigenvalues, "
      f"{theorems['open']} open gaps remaining.")
    p(sep)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="RH certificate ledger dashboard.")
    parser.add_argument("--source-dir", default="source", help="Source directory with .json files")
    parser.add_argument("--output", default=None, help="Write report to file (default: stdout)")
    parser.add_argument("--json-out", default=None, help="Write structured JSON summary to file")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    source_dir = (script_dir / args.source_dir).resolve()
    if not source_dir.is_dir():
        sys.exit(f"Source directory not found: {source_dir}")

    # Load all JSON files
    groups = {
        "theorem": [],
        "certificate": [],
        "scan": [],
        "other": [],
        "load_error": [],
    }

    total = 0
    for f in sorted(source_dir.glob("*.json")):
        total += 1
        try:
            data = json.loads(f.read_text(encoding="utf-8", errors="ignore"))
            cat = classify(f.name, data)
            groups[cat].append((f.name, data))
        except Exception as exc:
            groups["load_error"].append((f.name, {"_load_error": str(exc)}))

    print(f"Loaded {total} JSON files from {source_dir}", file=sys.stderr)

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        print_report(source_dir, groups, out)
    finally:
        if args.output:
            out.close()
            print(f"Report written: {args.output}", file=sys.stderr)

    if args.json_out:
        theorems = analyze_theorems(groups["theorem"])
        certs = analyze_certificates(groups["certificate"])
        scans, neg_any = analyze_scans(groups["scan"])
        summary = {
            "source_dir": str(source_dir),
            "totals": {k: len(v) for k, v in groups.items()},
            "theorems": theorems,
            "certificates": certs,
            "scans": scans,
            "neg_eigenvalue_scans": neg_any,
        }
        Path(args.json_out).write_text(json.dumps(summary, indent=2))
        print(f"JSON summary written: {args.json_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
