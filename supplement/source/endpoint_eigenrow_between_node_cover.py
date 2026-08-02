#!/usr/bin/env python3
r"""Parallel finite cover for exact between-node eigenrow disks."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import subprocess
import sys
from pathlib import Path


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": bool(closed),
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def run_one(command: list[str], output: str) -> dict:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return {
        "output": output,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def coverage(disks: list[dict], left: float, right: float) -> tuple[bool, float]:
    intervals = sorted(
        (float(row["realCell"][0]), float(row["realCell"][1]))
        for row in disks
        if row.get("radiiPolynomialClosed")
    )
    frontier = left
    for a, b in intervals:
        if a > frontier + 1e-15:
            return False, frontier
        frontier = max(frontier, b)
        if frontier >= right - 1e-15:
            return True, frontier
    return False, frontier


def public_disk(row: dict) -> dict:
    """Remove execution-cache paths from the theorem-level aggregate."""
    clean = json.loads(json.dumps(row))
    clean.pop("candidateSource", None)
    clean.pop("orientationCachePath", None)
    quadrature = clean.get("quadrature")
    if isinstance(quadrature, dict):
        quadrature.pop("cachePath", None)
    return clean


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", default="0.15125")
    parser.add_argument("--right", default="0.191808480511711849413426011005245")
    parser.add_argument("--spacing", type=float, default=0.0064)
    parser.add_argument("--degree", type=int, default=32)
    parser.add_argument("--inner-radius", default="0.006")
    parser.add_argument("--real-half-width", default="0.0033")
    parser.add_argument("--outer-radius", default="0.1")
    parser.add_argument("--quadrature-order", type=int, default=352)
    parser.add_argument("--bits", type=int, default=320)
    parser.add_argument("--majorant-step", default="0.02")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--run-dir", default=".endpoint_between_node_cover_runs")
    parser.add_argument(
        "--json-out", default="endpoint_eigenrow_between_node_crossing_cover.json"
    )
    args = parser.parse_args()

    left = float(args.left)
    right = float(args.right)
    count = max(1, math.ceil((right - left) / args.spacing))
    centers = [left + (right - left) * index / count for index in range(count + 1)]
    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    script = str(Path(__file__).with_name("endpoint_eigenrow_between_node_cauchy_certificate.py"))
    jobs = []
    half = float(args.real_half_width)
    for index, center in enumerate(centers):
        output = run_dir / f"disk_{index:03d}.json"
        command = [
            sys.executable,
            script,
            "--degree",
            str(args.degree),
            "--center",
            format(center, ".17g"),
            "--left",
            format(max(left, center - half), ".17g"),
            "--right",
            format(min(right, center + half), ".17g"),
            "--inner-radius",
            args.inner_radius,
            "--real-half-width",
            args.real_half_width,
            "--outer-radius",
            args.outer_radius,
            "--quadrature-order",
            str(args.quadrature_order),
            "--bits",
            str(args.bits),
            "--majorant-step",
            args.majorant_step,
            "--json-out",
            str(output),
        ]
        jobs.append((command, str(output)))

    print(
        f"between-node cover: {len(jobs)} disks, workers={args.workers}, "
        f"interval=[{args.left},{args.right}]",
        flush=True,
    )
    if args.reuse_existing:
        results = [
            {
                "output": output,
                "returncode": 0 if Path(output).exists() else 1,
                "stdout": "reused existing child certificate",
                "stderr": "" if Path(output).exists() else "child certificate is missing",
            }
            for _command, output in jobs
        ]
    else:
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(run_one, command, output) for command, output in jobs]
            for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
                row = future.result()
                results.append(row)
                print(
                    f"  completed {index}/{len(jobs)} rc={row['returncode']} output={row['output']}",
                    flush=True,
                )

    disks = []
    failures = []
    for row in results:
        path = Path(row["output"])
        if row["returncode"] != 0 or not path.exists():
            failures.append(row)
            continue
        child = json.loads(path.read_text(encoding="utf-8"))
        disks.extend(public_disk(disk) for disk in child.get("disks", []))
    covered, frontier = coverage(disks, left, right)
    local_closed = len(disks) == len(jobs) and all(
        row.get("radiiPolynomialClosed") for row in disks
    )
    theorem_closed = not failures and local_closed and covered
    data = {
        "schemaVersion": 1,
        "theoremName": "finite overlapping between-node eigenrow disk cover",
        "proofClass": "interval/ball certificate",
        "statement": (
            "The exact normalized negative endpoint eigenrow and its derivatives through "
            "order eight are enclosed on the requested interval by the union of the "
            "reported analytic Newton/Cauchy disks."
        ),
        "domain": {
            "parameterInterval": [left, right],
            "matrixFamily": "exact confluent endpoint integral J(s)",
        },
        "dependencies": [
            "endpoint_eigenrow_between_node_cauchy_certificate.py",
            "endpoint_confluent_segment_bernstein_certificate.json",
        ],
        "settings": {
            "degree": args.degree,
            "spacing": args.spacing,
            "innerRadius": args.inner_radius,
            "realHalfWidth": args.real_half_width,
            "outerRadius": args.outer_radius,
            "quadratureOrder": args.quadrature_order,
            "bits": args.bits,
            "majorantStep": args.majorant_step,
            "workers": args.workers,
        },
        "childCertificateCount": len(results),
        "failures": failures,
        "disks": sorted(disks, key=lambda row: row["center"]),
        "coverageFrontier": frontier,
        "allLocalDisksClosed": local_closed,
        "requestedIntervalCovered": covered,
        "betweenNodeEigenrowDerivativeCoverClosed": theorem_closed,
        "betweenNodeEigenrowDerivativeCoverStatus": status(
            "exact finite between-node eigenrow derivative cover",
            theorem_closed,
            (
                "Every local radii polynomial closes and the closed cells cover the interval."
                if theorem_closed
                else "At least one local disk failed or the certified cells leave a coverage gap."
            ),
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"cover closed: {theorem_closed}; frontier={frontier}")
    print(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
