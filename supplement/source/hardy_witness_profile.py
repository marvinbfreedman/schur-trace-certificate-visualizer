#!/usr/bin/env python3
"""Profile the Riesz witness w=A^+Bu for endpoint Hardy structure.

This reads douglas_top_vector.json and computes simple concentration metrics
for the top Douglas vector:

  - L2 mass fraction in right endpoint layers for u, source, and witness
  - L2 mass fraction in left endpoint layers for u, source, and witness
  - sup-norm locations
  - sign-change counts

These are not proof quantities by themselves, but they say what a Hardy/trace
estimate should control.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def trapz(points, key, lo=None, hi=None):
    total = 0.0
    for a, b in zip(points[:-1], points[1:]):
        mid = 0.5 * (a["s"] + b["s"])
        if lo is not None and mid < lo:
            continue
        if hi is not None and mid > hi:
            continue
        dx = b["s"] - a["s"]
        total += 0.5 * dx * (a[key] * a[key] + b[key] * b[key])
    return total


def max_abs(points, key):
    idx = max(range(len(points)), key=lambda i: abs(points[i][key]))
    return points[idx]["s"], points[idx][key]


def sign_changes(values):
    count = 0
    prev = 0
    for value in values:
        sign = 1 if value > 0 else -1 if value < 0 else 0
        if sign and prev and sign != prev:
            count += 1
        if sign:
            prev = sign
    return count


def fmt(value):
    if abs(value) >= 1e-3 and abs(value) < 1e4:
        return f"{value:.6g}"
    return f"{value:.3e}"


def layer_rows(points, key, length, layers):
    total = trapz(points, key)
    out = []
    for width in layers:
        left = trapz(points, key, 0.0, width)
        right = trapz(points, key, length - width, length)
        out.append(
            {
                "width": width,
                "left": left / total if total else 0.0,
                "right": right / total if total else 0.0,
            }
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="douglas_top_vector.json")
    parser.add_argument("--json-out", default="hardy_witness_profile.json")
    parser.add_argument("--layers", default="0.05 0.1 0.2 0.4 0.8")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    points = data["grid"]
    length = data["L"]
    layers = [float(piece) for piece in args.layers.replace(",", " ").split()]
    profile = {
        "source": args.input,
        "model": data["model"],
        "basis": data["basis"],
        "constraints": data["constraints"],
        "gamma2Top": data["gamma2Top"],
        "gammaTop": data["gammaTop"],
        "uPeak": dict(zip(("s", "value"), max_abs(points, "u"))),
        "sourcePeak": dict(zip(("s", "value"), max_abs(points, "source"))),
        "witnessPeak": dict(zip(("s", "value"), max_abs(points, "nWitness"))),
        "uSignChanges": sign_changes([p["u"] for p in points]),
        "sourceSignChanges": sign_changes([p["source"] for p in points]),
        "witnessSignChanges": sign_changes([p["nWitness"] for p in points]),
        "uLayers": layer_rows(points, "u", length, layers),
        "sourceLayers": layer_rows(points, "source", length, layers),
        "witnessLayers": layer_rows(points, "nWitness", length, layers),
    }
    Path(args.json_out).write_text(json.dumps(profile, indent=2), encoding="utf-8")

    print(
        f"Hardy witness profile model={profile['model']} basis={profile['basis']} "
        f"gamma2={fmt(profile['gamma2Top'])}"
    )
    print(
        f"  u_peak s={fmt(profile['uPeak']['s'])} value={fmt(profile['uPeak']['value'])} "
        f"sign_changes={profile['uSignChanges']}"
    )
    print(
        f"  source_peak s={fmt(profile['sourcePeak']['s'])} "
        f"value={fmt(profile['sourcePeak']['value'])} "
        f"sign_changes={profile['sourceSignChanges']}"
    )
    print(
        f"  witness_peak s={fmt(profile['witnessPeak']['s'])} "
        f"value={fmt(profile['witnessPeak']['value'])} "
        f"sign_changes={profile['witnessSignChanges']}"
    )
    print("  layer mass fractions:")
    print("    width     u_L      u_R      h_L      h_R      w_L      w_R")
    for u_row, h_row, w_row in zip(
        profile["uLayers"], profile["sourceLayers"], profile["witnessLayers"]
    ):
        print(
            f"    {u_row['width']:5.2f} "
            f"{u_row['left']:8.4f} {u_row['right']:8.4f} "
            f"{h_row['left']:8.4f} {h_row['right']:8.4f} "
            f"{w_row['left']:8.4f} {w_row['right']:8.4f}"
        )
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
