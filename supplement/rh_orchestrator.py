#!/usr/bin/env python3
"""
rh_orchestrator.py — Run RH verification scripts in parallel and collect results.

Usage:
    python rh_orchestrator.py [options]

Options:
    --source-dir PATH     Source directory (default: ./source)
    --workers N           Parallel workers (default: cpu_count - 1)
    --timeout SECS        Per-script timeout in seconds (default: 90)
    --filter PATTERN      Only run scripts whose name contains PATTERN
    --skip-existing       Skip scripts that already have a paired .json output
    --entry-points-only   Only run scripts not imported as libraries (default: True)
    --all                 Run all scripts including library modules
    --dry-run             Print scripts that would run, then exit
    --output FILE         Write JSON report to FILE (default: rh_run_report.json)
    --text-report FILE    Write text summary to FILE (default: rh_run_report.txt)
"""

import argparse
import concurrent.futures
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path


def find_entry_points(source_dir: Path) -> list:
    """Return scripts with __main__ guard that are not imported as libraries."""
    py_files = list(source_dir.glob("*.py"))

    imported = set()
    for f in py_files:
        try:
            text = f.read_text(encoding="utf-8", errors="ignore")
            for m in re.finditer(r'^(?:from|import)\s+([a-zA-Z_][a-zA-Z0-9_]*)', text, re.MULTILINE):
                imported.add(m.group(1))
        except Exception:
            pass

    entry_points = []
    for f in sorted(py_files):
        if f.stem in imported:
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="ignore")
            if "if __name__" in text and "__main__" in text:
                entry_points.append(f)
        except Exception:
            pass

    return entry_points


def find_all_runnable(source_dir: Path) -> list:
    """Return all scripts with __main__ guard."""
    result = []
    for f in sorted(source_dir.glob("*.py")):
        try:
            text = f.read_text(encoding="utf-8", errors="ignore")
            if "if __name__" in text and "__main__" in text:
                result.append(f)
        except Exception:
            pass
    return result


def has_existing_json(script: Path) -> bool:
    json_path = script.with_suffix(".json")
    return json_path.exists()


def run_script(task):
    script_path, source_dir, timeout = task
    start = time.monotonic()
    env = {**os.environ, "PYTHONPATH": str(source_dir)}
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(source_dir),
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.monotonic() - start
        return {
            "script": script_path.name,
            "returncode": result.returncode,
            "elapsed_s": round(elapsed, 3),
            "stdout": result.stdout[-3000:].strip() if result.stdout else "",
            "stderr": result.stderr[-1000:].strip() if result.stderr else "",
            "status": "ok" if result.returncode == 0 else "fail",
        }
    except subprocess.TimeoutExpired:
        return {
            "script": script_path.name,
            "returncode": None,
            "elapsed_s": round(time.monotonic() - start, 3),
            "stdout": "",
            "stderr": f"TIMEOUT after {timeout}s",
            "status": "timeout",
        }
    except Exception as exc:
        return {
            "script": script_path.name,
            "returncode": None,
            "elapsed_s": round(time.monotonic() - start, 3),
            "stdout": "",
            "stderr": str(exc),
            "status": "error",
        }


def write_text_report(results, output_path: Path, elapsed_total: float):
    ok = [r for r in results if r["status"] == "ok"]
    fail = [r for r in results if r["status"] == "fail"]
    timeout = [r for r in results if r["status"] == "timeout"]
    skipped = [r for r in results if r["status"] == "skipped"]
    error = [r for r in results if r["status"] == "error"]

    lines = [
        "=" * 70,
        "RH VERIFICATION ORCHESTRATOR — RUN REPORT",
        "=" * 70,
        f"Total scripts:  {len(results)}",
        f"  Passed:       {len(ok)}",
        f"  Failed:       {len(fail)}",
        f"  Timeout:      {len(timeout)}",
        f"  Skipped:      {len(skipped)}",
        f"  Error:        {len(error)}",
        f"Wall-clock time: {elapsed_total:.1f}s",
        "",
    ]

    if fail:
        lines.append("FAILURES:")
        for r in fail:
            lines.append(f"  {r['script']} (rc={r['returncode']}, {r['elapsed_s']:.1f}s)")
            if r["stderr"]:
                for ln in r["stderr"].splitlines()[-5:]:
                    lines.append(f"    stderr: {ln}")
        lines.append("")

    if timeout:
        lines.append("TIMEOUTS:")
        for r in timeout:
            lines.append(f"  {r['script']} ({r['elapsed_s']:.1f}s)")
        lines.append("")

    if error:
        lines.append("LAUNCH ERRORS:")
        for r in error:
            lines.append(f"  {r['script']}: {r['stderr']}")
        lines.append("")

    lines.append("PASSED:")
    for r in sorted(ok, key=lambda x: x["elapsed_s"], reverse=True):
        lines.append(f"  {r['elapsed_s']:6.2f}s  {r['script']}")
    lines.append("")

    output_path.write_text("\n".join(lines))
    print(f"Text report written: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Run RH verification scripts in parallel.")
    parser.add_argument("--source-dir", default="source", help="Source directory")
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1))
    parser.add_argument("--timeout", type=float, default=90.0, help="Per-script timeout (seconds)")
    parser.add_argument("--filter", dest="pattern", default=None, help="Name substring filter")
    parser.add_argument("--skip-existing", action="store_true", help="Skip scripts with paired .json")
    parser.add_argument("--all", dest="run_all", action="store_true", help="Run all scripts including libraries")
    parser.add_argument("--dry-run", action="store_true", help="Print scripts and exit")
    parser.add_argument("--output", default="rh_run_report.json", help="JSON report output path")
    parser.add_argument("--text-report", default="rh_run_report.txt", help="Text report output path")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    source_dir = (script_dir / args.source_dir).resolve()
    if not source_dir.is_dir():
        sys.exit(f"Source directory not found: {source_dir}")

    print(f"Source directory: {source_dir}")
    print(f"Workers: {args.workers}  Timeout: {args.timeout}s")

    if args.run_all:
        candidates = find_all_runnable(source_dir)
    else:
        candidates = find_entry_points(source_dir)
    print(f"Entry-point scripts found: {len(candidates)}")

    if args.pattern:
        candidates = [f for f in candidates if args.pattern in f.name]
        print(f"After filter '{args.pattern}': {len(candidates)}")

    results = []
    to_run = []
    for script in candidates:
        if args.skip_existing and has_existing_json(script):
            results.append({
                "script": script.name,
                "returncode": 0,
                "elapsed_s": 0.0,
                "stdout": "",
                "stderr": "",
                "status": "skipped",
            })
        else:
            to_run.append(script)

    if args.dry_run:
        print(f"\nWould run {len(to_run)} scripts (skipping {len(results)}):")
        for s in to_run:
            print(f"  {s.name}")
        return

    print(f"Scripts to run: {len(to_run)}  (skipped: {len(results)})")
    print()

    tasks = [(s, source_dir, args.timeout) for s in to_run]
    completed = 0
    wall_start = time.monotonic()

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_script, t): t[0].name for t in tasks}
        for fut in concurrent.futures.as_completed(futures):
            r = fut.result()
            results.append(r)
            completed += 1
            mark = "OK" if r["status"] == "ok" else r["status"].upper()
            print(f"  [{completed:4d}/{len(to_run)}] {mark:7s}  {r['elapsed_s']:6.2f}s  {r['script']}")
            if r["status"] in ("fail", "error") and r["stderr"]:
                for ln in r["stderr"].splitlines()[-3:]:
                    print(f"             {ln}")

    elapsed_total = time.monotonic() - wall_start
    print()

    ok_count = sum(1 for r in results if r["status"] == "ok")
    fail_count = sum(1 for r in results if r["status"] in ("fail", "timeout", "error"))
    skip_count = sum(1 for r in results if r["status"] == "skipped")
    print(f"Done in {elapsed_total:.1f}s — {ok_count} passed, {fail_count} failed/timeout, {skip_count} skipped")

    report = {
        "summary": {
            "total": len(results),
            "ok": ok_count,
            "failed": fail_count,
            "skipped": skip_count,
            "elapsed_s": round(elapsed_total, 2),
            "workers": args.workers,
            "timeout_s": args.timeout,
        },
        "results": sorted(results, key=lambda r: r["script"]),
    }
    Path(args.output).write_text(json.dumps(report, indent=2))
    print(f"JSON report written: {args.output}")

    write_text_report(results, Path(args.text_report), elapsed_total)


if __name__ == "__main__":
    main()
