#!/usr/bin/env python3
r"""Scan the Riesz spectral split under Galerkin refinement.

The single-profile script diagonalizes A on ker(R) for one finite section.
This wrapper repeats that computation and records whether the Douglas energy
is carried by a stable low-dimensional spectral window or by moving Galerkin
artifacts.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from riesz_spectral_profile import compute  # noqa: E402


def parse_ints(text):
    return [int(piece) for piece in text.replace(",", " ").split()]


def fmt(x, digits=8):
    return mp.nstr(x, digits)


def make_case(base, basis):
    if base.constraint_rule == "half":
        constraints = max(2, basis // 2)
    elif base.constraint_rule == "basis":
        constraints = basis
    else:
        constraints = base.constraints
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=base.quad_base + base.quad_step * basis,
        laguerre=base.laguerre,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
        dps=base.dps,
        points=base.points,
        top=base.top,
    )


def summarize(profile):
    rows = profile["rows"]
    top_energy = profile["topEnergyModes"]
    top_witness = profile["topWitnessModes"]
    energy_low2 = sum(row["energyFrac"] for row in rows[: min(2, len(rows))])
    witness_low2 = sum(row["witnessFrac"] for row in rows[: min(2, len(rows))])
    energy_top3 = sum(row["energyFrac"] for row in top_energy[:3])
    energy_eff_rank = 1.0 / sum(row["energyFrac"] ** 2 for row in rows if row["energyFrac"])
    return {
        "basis": profile["basis"],
        "constraints": profile["constraints"],
        "gamma2Top": profile["gamma2Top"],
        "kerMin": profile["kerMin"],
        "energyLow2": energy_low2,
        "witnessLow2": witness_low2,
        "energyTop3": energy_top3,
        "energyEffectiveRank": energy_eff_rank,
        "topEnergyModes": [
            {
                "mode": row["mode"],
                "energyFrac": row["energyFrac"],
                "witnessFrac": row["witnessFrac"],
                "lambda": row["lambda"],
                "modePeakS": row["modePeakS"],
            }
            for row in top_energy[:3]
        ],
        "topWitnessModes": [
            {
                "mode": row["mode"],
                "energyFrac": row["energyFrac"],
                "witnessFrac": row["witnessFrac"],
                "lambda": row["lambda"],
                "modePeakS": row["modePeakS"],
            }
            for row in top_witness[:3]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="endpoint_b")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", default="10 12 14 16")
    parser.add_argument("--constraint-rule", choices=["half", "fixed", "basis"], default="half")
    parser.add_argument("--constraints", type=int, default=8)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--quad-base", type=int, default=6)
    parser.add_argument("--quad-step", type=int, default=1)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--points", type=int, default=121)
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--json-out", default="riesz_spectral_scan.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    summaries = []
    print(
        f"Riesz spectral scan model={args.model} kind={args.kind} "
        f"omega={args.omega} L={args.L}",
        flush=True,
    )
    print(
        "  basis cons gamma2       low2_E  low2_W  top3_E  eff_rank  top_E_modes",
        flush=True,
    )
    for basis in parse_ints(args.basis):
        case = make_case(args, basis)
        profile = compute(case)
        summary = summarize(profile)
        summaries.append(summary)
        modes = ",".join(str(row["mode"]) for row in summary["topEnergyModes"])
        print(
            f"  {summary['basis']:5d} {summary['constraints']:4d} "
            f"{fmt(summary['gamma2Top']):>10} "
            f"{summary['energyLow2']:7.3f} {summary['witnessLow2']:7.3f} "
            f"{summary['energyTop3']:7.3f} {summary['energyEffectiveRank']:8.3f} "
            f"{modes}",
            flush=True,
        )

    data = {
        "model": args.model,
        "kind": args.kind,
        "omega": float(mp.mpf(args.omega)),
        "L": float(mp.mpf(args.L)),
        "rows": summaries,
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
