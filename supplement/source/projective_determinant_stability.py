#!/usr/bin/env python3
"""Projective/column-normalized determinant stability constants.

The raw determinant perturbation bound compares ``C alpha`` with a raw
singular value.  The derivative-rank certificate, however, uses the matrix
whose columns have been normalized.  This script converts the raw diagnostic
into the scale-aware projective bound:

    ||hat M_N - hat M|| <= 2 sqrt(d) C alpha / d_min,

where d_min is the smallest derivative-response column norm.  Thus the
normalized determinant remains nonzero if

    alpha < eta d_min / (2 sqrt(d) C),

where eta is the normalized singular margin.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="determinant_gap_bound_diagnostic.json")
    parser.add_argument("--json-out", default="projective_determinant_stability.json")
    args = parser.parse_args()

    source = json.loads(Path(args.input).read_text(encoding="utf-8"))
    rows = []
    for row in source["rows"]:
        dim = int(row["activeDim"])
        eta = float(row["normalizedMinSingularValue"])
        c_high = float(row["COnHighBlock"])
        col_norms = [float(x) for x in row["columnNorms"]]
        d_min = min(col_norms) if col_norms else 0.0
        d_max = max(col_norms) if col_norms else 0.0
        scale_ratio = d_max / d_min if d_min else math.inf
        projective_alpha = eta * d_min / (2 * math.sqrt(dim) * c_high) if c_high and d_min else 0.0
        half_margin_alpha = projective_alpha / 2
        rows.append(
            {
                "basis": int(row["basis"]),
                "activeDim": dim,
                "normalizedMargin": eta,
                "minColumnNorm": d_min,
                "maxColumnNorm": d_max,
                "columnScaleRatio": scale_ratio,
                "COnHighBlock": c_high,
                "projectiveAllowedAlpha": projective_alpha,
                "halfMarginProjectiveAllowedAlpha": half_margin_alpha,
                "rawAllowedAlpha": float(row["allowedAlphaRaw"]),
            }
        )

    min_projective = min((row["projectiveAllowedAlpha"] for row in rows), default=0.0)
    max_scale_ratio = max((row["columnScaleRatio"] for row in rows), default=0.0)
    data = {
        "source": args.input,
        "point": source.get("point"),
        "orders": source.get("orders"),
        "bases": source.get("bases"),
        "minProjectiveAllowedAlpha": min_projective,
        "maxColumnScaleRatio": max_scale_ratio,
        "rows": rows,
        "interpretation": (
            "Column normalization is legitimate only with explicit column-scale "
            "stability.  The projective allowed alpha is the sufficient "
            "subspace/eigenvector perturbation size preserving the normalized "
            "singular margin by the elementary normalization Lipschitz bound."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Projective determinant stability")
    print("  basis dim eta d_min/C scale_ratio allowed_alpha")
    for row in rows:
        dmin_over_c = row["minColumnNorm"] / row["COnHighBlock"] if row["COnHighBlock"] else 0.0
        print(
            "  {basis:5d} {activeDim:3d} {eta:.6g} {dmin_over_c:.6g} "
            "{columnScaleRatio:.6g} {projectiveAllowedAlpha:.6g}".format(
                eta=row["normalizedMargin"],
                dmin_over_c=dmin_over_c,
                **row,
            )
        )
    print(f"  min projective allowed alpha: {min_projective:.6g}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
