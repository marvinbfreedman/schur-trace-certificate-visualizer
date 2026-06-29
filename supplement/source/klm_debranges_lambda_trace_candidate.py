#!/usr/bin/env python3
r"""Apply Lambda_a to primitive de Branges/Weyl evaluation candidates.

For a primitive candidate f_z(s), the finite sampled trace is

    y_z(a_j)=Lambda_{a_j}(f_z)
      = sum_k e_k(a_j) f_z^(k)(a_j)/k!.

The Volterra Green certificates use row-space coordinates c_z rather than raw
sample values.  If U spans the row-space of R^*R, then

    R U c_z ~= y_z.

This script computes c_z for several exponential primitive candidates and
checks the resulting Volterra/KLM branch Grams against the canonical Hardy
branch Grams.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    block_residuals,
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    herm,
    local_args,
    mapped_joint,
    rel_frob,
    to_float,
)
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from quotient_factorization_mp import sampled_centers  # noqa: E402


def lambda_row(center: mp.mpf, args: argparse.Namespace):
    vals, neg, vec = defect_at(
        center,
        args.jet_order,
        mp.mpf(args.endpoint_rmax),
        args.endpoint_order,
        mp.mpf(args.endpoint_tol),
    )
    if vec[0] > 0:
        vec = [-x for x in vec]
    return vec


def candidate_alpha(z: complex, kind: str, omega: float):
    if kind == "exp_iz":
        return 1j * z, 1.0 + 0j
    if kind == "exp_minus_iz":
        return -1j * z, 1.0 + 0j
    if kind == "plus_weighted":
        return 1j * z, None
    if kind == "minus_weighted":
        return 1j * z, None
    if kind == "omega_shift_plus":
        return 1j * z - omega, 1.0 + 0j
    if kind == "omega_shift_minus":
        return 1j * z + omega, 1.0 + 0j
    raise ValueError(kind)


def lambda_of_exp(center: mp.mpf, row, alpha: complex, amp: complex) -> mp.mpc:
    a = mp.mpc(alpha)
    total = mp.mpc(0)
    for k, coeff in enumerate(row):
        total += coeff * (a ** k) * mp.e ** (a * center) / mp.factorial(k)
    return mp.mpc(amp) * total


def solve_trace_coordinates(model: dict, y: mp.matrix) -> mp.matrix:
    A = model["R"] * model["U"]
    lhs = ctranspose(A) * A
    rhs = ctranspose(A) * y
    reg = mp.mpf("1e-40")
    for i in range(min(lhs.rows, lhs.cols)):
        lhs[i, i] += reg
    return mp.lu_solve(lhs, rhs)


def candidate_coordinates(args: argparse.Namespace, model: dict, xi: XiTransform, kind: str):
    centers = sampled_centers(
        mp.mpf(args.constraint_min),
        mp.mpf(args.constraint_max),
        args.constraints,
    )
    rows = [lambda_row(center, args) for center in centers]
    nodes = [
        0.0 + 0.35j,
        0.45 + 0.40j,
        0.95 + 0.55j,
        1.55 + 0.75j,
        2.25 + 0.95j,
    ]
    X = mp.matrix(model["traceRank"], len(nodes))
    omega = float(args.omega)
    for j, z in enumerate(nodes):
        alpha, amp = candidate_alpha(z, kind, omega)
        if kind == "plus_weighted":
            amp = xi(z + 1j * omega)
        elif kind == "minus_weighted":
            amp = xi(z - 1j * omega)
        y = mp.matrix(len(centers), 1)
        for i, center in enumerate(centers):
            y[i] = lambda_of_exp(center, rows[i], alpha, amp)
        c = solve_trace_coordinates(model, y)
        for i in range(c.rows):
            X[i, j] = c[i]
    return X


def fit_global_scale(mapped: mp.matrix, target: mp.matrix) -> tuple[mp.mpf, mp.mpf]:
    num = mp.mpf("0")
    den = mp.mpf("0")
    for i in range(mapped.rows):
        for j in range(mapped.cols):
            num += mp.re(mp.conj(mapped[i, j]) * target[i, j])
            den += abs(mapped[i, j]) ** 2
    scale = max(mp.mpf("0"), num / den) if den else mp.mpf("0")
    return scale, rel_frob(scale * mapped, target)


def evaluate_candidate(args, model, target, X):
    mapped = mapped_joint(model, X)
    n = target["P"].rows
    raw = {k: to_float(v) for k, v in block_residuals(mapped, target["joint"], n).items()}
    scale, scaled_joint = fit_global_scale(mapped, target["joint"])
    raw["globalScale"] = to_float(scale)
    raw["scaledJointRelative"] = to_float(scaled_joint)
    raw["traceCoordinateNormFrobenius"] = to_float(frob_norm(X))
    raw["mappedJointNorm"] = to_float(frob_norm(mapped))
    return raw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=5)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-24")
    parser.add_argument("--psd-tol", default="1e-24")
    parser.add_argument("--s-order", type=int, default=20)
    parser.add_argument("--u-order", type=int, default=40)
    parser.add_argument("--u-max", default="10")
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=800)
    parser.add_argument("--hardy-rmax", type=float, default=40.0)
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--json-out", default="klm_debranges_lambda_trace_candidate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    model = finite_trace_model(args)
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    full = hardy_targets(xi, float(args.omega), None)
    trunc = hardy_targets(xi, float(args.omega), args.hardy_rmax)
    candidate_kinds = [
        "exp_iz",
        "exp_minus_iz",
        "omega_shift_plus",
        "omega_shift_minus",
        "plus_weighted",
        "minus_weighted",
    ]
    results = []
    for kind in candidate_kinds:
        X = candidate_coordinates(args, model, xi, kind)
        results.append(
            {
                "kind": kind,
                "full": evaluate_candidate(args, model, full, X),
                "truncated": evaluate_candidate(args, model, trunc, X),
            }
        )
    best = min(results, key=lambda row: row["full"]["scaledJointRelative"])
    data = {
        "theoremName": "Lambda trace candidate for de Branges evaluation vector",
        "basis": args.basis,
        "constraints": args.constraints,
        "traceRank": model["traceRank"],
        "traceNullity": model["traceNullity"],
        "omega": float(args.omega),
        "candidateFormula": (
            "y_z(a)=Lambda_a(f_z), c_z solves R U c_z=y_z, then x_z=c_z "
            "is tested in the Volterra/KLM branch Gram model."
        ),
        "results": results,
        "bestByScaledFullJointResidual": best,
        "exactTraceMapConstructed": bool(best["full"]["scaledJointRelative"] < 1e-6),
        "diagnosis": (
            "None of the elementary exponential primitive candidates produces "
            "the joint Hardy branch Grams.  The primitive/evaluation vector "
            "must include the missing Weyl/Volterra integral transform, not "
            "just the bare exponential jet."
            if best["full"]["scaledJointRelative"] >= 1e-6
            else "An elementary Lambda trace candidate matches the joint branch Grams."
        ),
        "nextProofTarget": (
            "Derive the primitive evaluation vector by inverting the Volterra "
            "feature map: solve G_+(f_z)=h_z^+ at the feature level, then apply "
            "Lambda_a to that f_z."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Lambda trace candidate for de Branges evaluation vector")
    print(f"  basis={args.basis} constraints={args.constraints} trace_rank={model['traceRank']}")
    for row in results:
        full = row["full"]
        print(
            f"  {row['kind']:16s} raw={full['jointRelative']:.6e} "
            f"scaled={full['scaledJointRelative']:.6e} scale={full['globalScale']:.3e}"
        )
    print(f"  best={best['kind']} scaled={best['full']['scaledJointRelative']:.6e}")
    print(f"  exact trace map constructed: {data['exactTraceMapConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
