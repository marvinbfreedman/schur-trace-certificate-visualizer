#!/usr/bin/env python3
r"""Refinement scan for the source-window regularizer certificate.

The single source-window certificate is only useful if its tiny high/full
fraction is stable under two numerical refinements:

  1. the closed-trace Galerkin basis;
  2. the source-window quadrature count.

This driver reruns source_window_regularizer_scan.py over a small grid of
`basis x source_count` cases and aggregates the output.  The default basis
rule follows the earlier scans:

    constraints = floor(0.625 * basis).

The cutoff remains M=6, so the reported fractions are for the high block
H_6^\perp inside the closed trace space.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    ax = abs(x)
    if 1e-3 <= ax < 1e4:
        return f"{x:.6g}"
    return f"{x:.6e}"


def constraints_for(args: argparse.Namespace, basis: int) -> int:
    if args.constraint_rule == "fixed":
        return args.constraints
    if args.constraint_rule == "ratio-floor":
        return max(1, min(basis - 1, math.floor(args.constraint_ratio * basis)))
    if args.constraint_rule == "ratio-round":
        return max(1, min(basis - 1, round(args.constraint_ratio * basis)))
    if args.constraint_rule == "offset":
        return max(1, min(basis - 1, basis - args.constraint_offset))
    raise ValueError(args.constraint_rule)


def run_child(args: argparse.Namespace, basis: int, source_count: int, out_path: Path) -> dict:
    constraints = constraints_for(args, basis)
    cmd = [
        sys.executable,
        "source_window_regularizer_scan.py",
        "--basis",
        str(basis),
        "--constraints",
        str(constraints),
        "--source-count",
        str(source_count),
        "--cutoff",
        str(args.cutoff),
        "--dps",
        str(args.dps),
        "--source-min",
        args.source_min,
        "--source-max",
        args.source_max,
        "--omega",
        args.omega,
        "--L",
        args.L,
        "--quad",
        str(args.quad),
        "--laguerre",
        str(args.laguerre),
        "--constraint-min",
        args.constraint_min,
        "--constraint-max",
        args.constraint_max,
        "--kernel-order",
        str(args.kernel_order),
        "--kernel-rmax",
        args.kernel_rmax,
        "--matrix-order",
        str(args.matrix_order),
        "--matrix-rmax",
        args.matrix_rmax,
        "--endpoint-order",
        str(args.endpoint_order),
        "--endpoint-rmax",
        args.endpoint_rmax,
        "--json-out",
        str(out_path),
    ]
    completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
    data = json.loads(out_path.read_text(encoding="utf-8"))
    data["stdout"] = completed.stdout
    return data


def summarize(data: dict) -> dict:
    return {
        "basis": data["basis"],
        "constraints": data["constraints"],
        "rank": data["rank"],
        "nullity": data["nullity"],
        "positiveModes": data["positiveModes"],
        "cutoff": data["cutoff"],
        "sourceCount": data["sourceCount"],
        "lambdaMinA": data["lambdaMinA"],
        "lambdaMaxA": data["lambdaMaxA"],
        "integratedHighAbsorb": data["integratedHighAbsorb"],
        "integratedFullAbsorb": data["integratedFullAbsorb"],
        "integratedHighFullFrac": data["integratedHighFullFrac"],
        "worstU": data["worstSingle"]["u"],
        "worstSingleHighAbsorb": data["worstSingle"]["highAbsorb"],
        "worstSingleFullAbsorb": data["worstSingle"]["fullAbsorb"],
        "worstSingleHighFullFrac": data["worstSingle"]["highFullFrac"],
        "emptyHighBlock": data["positiveModes"] <= data["cutoff"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", default="18 20 22")
    parser.add_argument("--source-counts", default="5 9 13")
    parser.add_argument("--constraint-rule", choices=["fixed", "ratio-floor", "ratio-round", "offset"], default="ratio-floor")
    parser.add_argument("--constraint-ratio", type=float, default=0.625)
    parser.add_argument("--constraint-offset", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=12)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--json-out", default="source_window_refinement_scan.json")
    args = parser.parse_args()

    bases = parse_ints(args.bases)
    counts = parse_ints(args.source_counts)
    rows = []
    full_cases = []

    print(
        "Source-window refinement scan "
        f"bases={bases} source_counts={counts} rule={args.constraint_rule}",
        flush=True,
    )
    print("  basis cons pos src integrated_frac worst_frac worst_u", flush=True)

    with tempfile.TemporaryDirectory(prefix="source-window-refinement-") as tmp:
        tmp_dir = Path(tmp)
        for basis in bases:
            for count in counts:
                out_path = tmp_dir / f"b{basis}_q{count}.json"
                data = run_child(args, basis, count, out_path)
                row = summarize(data)
                rows.append(row)
                full_cases.append(
                    {
                        key: value
                        for key, value in data.items()
                        if key not in {"stdout", "perSource"}
                    }
                )
                print(
                    f"  {basis:5d} {row['constraints']:4d} {row['positiveModes']:3d} "
                    f"{count:3d} {fmt(row['integratedHighFullFrac']):>15} "
                    f"{fmt(row['worstSingleHighFullFrac']):>11} "
                    f"{row['worstU']:.6f}",
                    flush=True,
                )

    frac_values = [row["integratedHighFullFrac"] for row in rows if not row["emptyHighBlock"]]
    worst_values = [row["worstSingleHighFullFrac"] for row in rows if not row["emptyHighBlock"]]
    summary = {
        "integratedMin": min(frac_values) if frac_values else None,
        "integratedMax": max(frac_values) if frac_values else None,
        "integratedSpread": (max(frac_values) - min(frac_values)) if frac_values else None,
        "worstSingleMin": min(worst_values) if worst_values else None,
        "worstSingleMax": max(worst_values) if worst_values else None,
        "worstSingleSpread": (max(worst_values) - min(worst_values)) if worst_values else None,
    }
    result = {
        "omega": float(args.omega),
        "L": float(args.L),
        "sourceMin": float(args.source_min),
        "sourceMax": float(args.source_max),
        "bases": bases,
        "sourceCounts": counts,
        "constraintRule": args.constraint_rule,
        "constraintRatio": args.constraint_ratio,
        "cutoff": args.cutoff,
        "summary": summary,
        "rows": rows,
        "cases": full_cases,
        "interpretation": (
            "Refines the source-window regularizer over Galerkin basis and "
            "Gauss source count. Stability of the high/full fraction supports "
            "a continuum source-window Hardy/Green estimate rather than a "
            "finite quadrature artifact."
        ),
    }
    Path(args.json_out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        "  summary integrated range "
        f"[{fmt(summary['integratedMin'])}, {fmt(summary['integratedMax'])}]",
        flush=True,
    )
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
