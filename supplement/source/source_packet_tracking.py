#!/usr/bin/env python3
"""Track the source-envelope packet behind the high-frequency tail estimate.

The source-envelope certificate orders modes by the positive operator A:

    A e_j = lambda_j e_j,
    s_j = lambda_j^{-1} ||b_j||^2 / C_full.

As the Galerkin section grows, new very small lambda_j can be inserted at the
front.  A fixed mode cutoff can therefore look unstable even when the source
packet is stable in spectral variables.  This script records both pictures:
mode-index packet summaries and tails above fixed spectral thresholds.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_THRESHOLDS = "1e-6 1e-5 1e-4 1e-3 1e-2 1e-1"


def parse_floats(text: str) -> list[float]:
    return [float(piece) for piece in text.replace(",", " ").split()]


def tail_at(row: dict, start: int) -> dict:
    tails = row["tails"]
    if start < len(tails):
        return tails[start]
    return tails[-1]


def first_tail_cutoff(row: dict, key: str, threshold: float) -> int | None:
    for tail in row["tails"]:
        if tail[key] <= threshold:
            return tail["start"]
    return None


def lambda_at_start(row: dict, start: int | None) -> float | None:
    if start is None:
        return None
    envelope = row["envelope"]
    if start >= len(envelope):
        return None
    return envelope[start]["lambda"]


def spectral_tail(row: dict, threshold: float) -> dict:
    envelope = row["envelope"]
    start = len(envelope)
    for item in envelope:
        if item["lambda"] >= threshold:
            start = item["mode"]
            break
    tail = tail_at(row, start)
    return {
        "lambdaThreshold": threshold,
        "start": start,
        "operatorFrac": tail["operatorFrac"],
        "scalarFrac": tail["scalarFrac"],
    }


def summarize_row(row: dict, significant_threshold: float, spectral_thresholds: list[float]) -> dict:
    envelope = row["envelope"]
    shares = [item["scalarFrac"] for item in envelope]
    total = sum(shares)
    peak = max(envelope, key=lambda item: item["scalarFrac"])

    if total > 0:
        center_mode = sum(item["mode"] * item["scalarFrac"] for item in envelope) / total
        width_mode = math.sqrt(
            sum(((item["mode"] - center_mode) ** 2) * item["scalarFrac"] for item in envelope) / total
        )
        center_log_lambda = sum(math.log(item["lambda"]) * item["scalarFrac"] for item in envelope) / total
        width_log_lambda = math.sqrt(
            sum(
                ((math.log(item["lambda"]) - center_log_lambda) ** 2) * item["scalarFrac"]
                for item in envelope
            )
            / total
        )
        lambda_geometric_center = math.exp(center_log_lambda)
    else:
        center_mode = 0.0
        width_mode = 0.0
        width_log_lambda = 0.0
        lambda_geometric_center = 0.0

    significant = [item for item in envelope if item["scalarFrac"] >= significant_threshold]
    first_sig = significant[0]["mode"] if significant else None
    last_sig = significant[-1]["mode"] if significant else None
    after_sig_start = (last_sig + 1) if last_sig is not None else 0
    after_sig_tail = tail_at(row, after_sig_start)

    scalar05_start = first_tail_cutoff(row, "scalarFrac", 0.05)
    operator05_start = first_tail_cutoff(row, "operatorFrac", 0.05)

    return {
        "basis": row["basis"],
        "constraints": row["constraints"],
        "rank": row["rank"],
        "nullity": row["nullity"],
        "positiveModes": row["positiveModes"],
        "gamma2": row["gamma2"],
        "totalScalarFrac": total,
        "peakMode": peak["mode"],
        "peakLambda": peak["lambda"],
        "peakShare": peak["scalarFrac"],
        "centerMode": center_mode,
        "widthMode": width_mode,
        "lambdaGeometricCenter": lambda_geometric_center,
        "widthLogLambda": width_log_lambda,
        "significantThreshold": significant_threshold,
        "firstSignificantMode": first_sig,
        "lastSignificantMode": last_sig,
        "significantModes": [item["mode"] for item in significant],
        "tailAfterSignificantStart": after_sig_start,
        "tailAfterSignificantOperatorFrac": after_sig_tail["operatorFrac"],
        "tailAfterSignificantScalarFrac": after_sig_tail["scalarFrac"],
        "tailCutoffScalar05": scalar05_start,
        "tailCutoffOperator05": operator05_start,
        "lambdaCutoffScalar05": lambda_at_start(row, scalar05_start),
        "lambdaCutoffOperator05": lambda_at_start(row, operator05_start),
        "spectralTails": [spectral_tail(row, threshold) for threshold in spectral_thresholds],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="source_envelope_tail_certificate.json")
    parser.add_argument("--json-out", default="source_packet_tracking.json")
    parser.add_argument("--significant-threshold", type=float, default=0.05)
    parser.add_argument("--spectral-thresholds", default=DEFAULT_THRESHOLDS)
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    thresholds = parse_floats(args.spectral_thresholds)
    rows = [
        summarize_row(row, args.significant_threshold, thresholds)
        for row in data["rows"]
    ]

    out = {
        "source": args.input,
        "model": data.get("model"),
        "kind": data.get("kind"),
        "omega": data.get("omega"),
        "L": data.get("L"),
        "significantThreshold": args.significant_threshold,
        "spectralThresholds": thresholds,
        "rows": rows,
    }
    Path(args.json_out).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    print("Source packet tracking")
    print("  basis peak(mode,lambda,share)  sig-modes tail-after-sig  tail(lambda>=1e-2)")
    for row in rows:
        tail_1e2 = next(
            item for item in row["spectralTails"]
            if abs(item["lambdaThreshold"] - 1e-2) < 1e-15
        )
        print(
            f"  {row['basis']:5d} "
            f"{row['peakMode']:3d} {row['peakLambda']:.3e} {row['peakShare']:.3f} "
            f"{','.join(str(m) for m in row['significantModes']):>12} "
            f"{row['tailAfterSignificantScalarFrac']:.3f} "
            f"{tail_1e2['scalarFrac']:.3f}"
        )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
