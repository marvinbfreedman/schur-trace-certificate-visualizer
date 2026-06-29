#!/usr/bin/env python3
"""Summarize active trace-response determinant certificates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_rows(paths):
    rows = []
    for path in paths:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        for row in data.get("rows", []):
            out = dict(row)
            out["sourceFile"] = path
            rows.append(out)
    rows.sort(key=lambda row: (int(row["basis"]), int(row["traceCount"])))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inputs",
        nargs="+",
        default=[
            "trace_active_minor_scan_t9.json",
            "trace_active_minor_scan_t13.json",
            "trace_active_minor_certificate.json",
        ],
    )
    parser.add_argument("--json-out", default="trace_active_minor_summary.json")
    args = parser.parse_args()

    rows = load_rows(args.inputs)
    seen = set()
    unique = []
    for row in rows:
        key = (int(row["basis"]), int(row["traceCount"]))
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)

    summary_rows = []
    for row in unique:
        minor = row.get("bestNormalizedMinor") or {}
        frame_vals = row.get("weightedFrameEigenvalues") or []
        summary_rows.append(
            {
                "basis": int(row["basis"]),
                "traceCount": int(row["traceCount"]),
                "activeDim": int(row["activeDim"]),
                "frameMin": float(frame_vals[0]) if frame_vals else 0.0,
                "bestAbsMinor": float(minor.get("absDeterminant", 0.0)),
                "traceCenters": minor.get("traceCenters", []),
                "sourceFile": row["sourceFile"],
            }
        )

    min_minor = min(row["bestAbsMinor"] for row in summary_rows) if summary_rows else 0.0
    data = {
        "inputs": args.inputs,
        "rowCount": len(summary_rows),
        "minBestAbsMinor": min_minor,
        "rows": summary_rows,
        "interpretation": (
            "Nonzero active trace-response minors prove finite sampled "
            "injectivity.  The continuum theorem should prove nonvanishing of "
            "the exterior trace determinant on the limiting active space."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Active trace minor summary")
    for row in summary_rows:
        print(
            "  basis {basis}, traces {traceCount}, dim {activeDim}: "
            "minor={bestAbsMinor:.6g}, centers={traceCenters}".format(**row)
        )
    print(f"  min best minor: {min_minor:.6g}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
