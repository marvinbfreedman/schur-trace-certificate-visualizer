#!/usr/bin/env python3
r"""Compose the high-block Hardy graph estimate with commuted kernel coercivity.

Inputs:

    ||E_u f||^2 <= H_{m,M}(u) ||f||_{W^{m,2}}^2,
    ||f||_{W^{m,2}}^2 <= G_{m,M} <S_m f,f>.

This script reports the product

    D_{m,M} = max_u H_{m,M}(u) G_{m,M},

which is the finite high-frequency domination constant

    ||E_u f||^2 <= D_{m,M} <S_m f,f>,     f in H_M cap ker R.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def fmt(x: float) -> str:
    return f"{x:.6g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hardy", default="lagrange_hardy_graph_certificate.json")
    parser.add_argument("--commuted", default="lagrange_commuted_kernel_energy.json")
    parser.add_argument("--json-out", default="lagrange_commuted_dominance_summary.json")
    args = parser.parse_args()

    hardy = json.loads(Path(args.hardy).read_text(encoding="utf-8"))
    commuted = json.loads(Path(args.commuted).read_text(encoding="utf-8"))

    hardy_max = {}
    hardy_arg = {}
    for row in hardy["rows"]:
        cutoff = int(row["cutoff"])
        t = float(row.get("t", math.nan))
        for item in row["graph"]:
            key = (cutoff, int(item["order"]))
            value = float(item["constant"])
            if key not in hardy_max or value > hardy_max[key]:
                hardy_max[key] = value
                hardy_arg[key] = t

    commuted_map = {
        (int(row["cutoff"]), int(row["order"])): row
        for row in commuted["rows"]
    }

    rows = []
    for key in sorted(set(hardy_max) & set(commuted_map)):
        cutoff, order = key
        hconst = hardy_max[key]
        gconst = float(commuted_map[key]["constant"])
        product = hconst * gconst
        rows.append(
            {
                "cutoff": cutoff,
                "order": order,
                "hardyMax": hconst,
                "hardyWorstT": hardy_arg[key],
                "commutedConstant": gconst,
                "product": product,
                "normProduct": math.sqrt(max(0.0, product)),
                "subunit": product < 1.0,
            }
        )

    print("Commuted high-block domination")
    print("  cutoff order   HardyMax   CommConst    Product    sqrt(Product)")
    for row in rows:
        print(
            f"  {row['cutoff']:6d} {row['order']:5d} "
            f"{fmt(row['hardyMax']):>10} {fmt(row['commutedConstant']):>11} "
            f"{fmt(row['product']):>10} {fmt(row['normProduct']):>13}"
        )

    data = {
        "basis": hardy["basis"],
        "constraints": hardy["constraints"],
        "omega": hardy["omega"],
        "L": hardy["L"],
        "rows": rows,
        "interpretation": (
            "Product of the high-block Hardy graph constant and the commuted "
            "Sturm elliptic constant. Product below one means the commuted "
            "energy dominates that residual row in the finite certificate."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
