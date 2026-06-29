#!/usr/bin/env python3
r"""Invert the finite Volterra plus-feature map G_+(f_z)=h_z^+.

This is the coefficient-level version of the trace-map construction.  We solve

    G_+ f_z ~= h_z^+

in the Galerkin basis, compute raw endpoint traces

    x_z = R f_z = Lambda(f_z),

then compare two induced branch Gram families with the Hardy branch Grams:

1. Direct feature Grams from the solved coefficients f_z.
2. Green-minimizer trace Grams obtained by converting x_z into trace-row
   coordinates and applying the completed Volterra trace model.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import default_z_nodes, hardy_feature_scalar  # noqa: E402
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    block2,
    block_residuals,
    branch_feature_matrices_mp,
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    herm,
    local_args,
    mapped_joint,
    rel_frob,
    stack_rows,
    sym,
    to_float,
)
from klm_debranges_lambda_trace_candidate import solve_trace_coordinates  # noqa: E402
from quotient_factorization_mp import shifted_legendre_polys, trace_matrix  # noqa: E402


def full_feature_model(args: argparse.Namespace):
    sigmas = [0] if mp.mpf(args.omega) == 0 else [-1, 1]
    branch_weight = mp.mpf("1") if mp.mpf(args.omega) == 0 else mp.mpf("0.5")
    plus_rows = None
    minus_rows = None
    row_us = []
    row_sigmas = []
    row_weights = []
    for sigma in sigmas:
        m, n, u_pts = branch_feature_matrices_mp(args, sigma)
        scale = mp.sqrt(branch_weight / 4)
        pblock = scale * (m + n)
        mblock = scale * (m - n)
        if plus_rows is None:
            plus_rows = pblock
            minus_rows = mblock
        else:
            plus_rows = stack_rows(plus_rows, pblock)
            minus_rows = stack_rows(minus_rows, mblock)
        row_us.extend(u_pts)
        row_sigmas.extend([sigma for _ in u_pts])
        row_weights.extend([branch_weight for _ in u_pts])
    polys = shifted_legendre_polys(args.basis, mp.mpf(args.L))
    _centers, R = trace_matrix(polys, local_args(args))
    P = sym(plus_rows.T * plus_rows)
    M = sym(minus_rows.T * minus_rows)
    C = plus_rows.T * minus_rows
    return {
        "plusRows": plus_rows,
        "minusRows": minus_rows,
        "rowUs": row_us,
        "rowSigmas": row_sigmas,
        "rowWeights": row_weights,
        "P": P,
        "M": M,
        "C": C,
        "R": R,
    }


def target_vector(args, xi: XiTransform, z: complex, model: dict, variant: str):
    amp = hardy_feature_scalar(xi, float(args.omega), z, "plus")
    out = mp.matrix(len(model["rowUs"]), 1)
    for i, u in enumerate(model["rowUs"]):
        sigma = model["rowSigmas"][i]
        base = amp * mp.e ** (1j * z * u)
        if variant == "same":
            value = base
        elif variant == "split":
            value = base / mp.sqrt(2) if sigma != 0 else base
        elif variant == "branch_unweighted":
            value = base / mp.sqrt(model["rowWeights"][i])
        elif variant == "undo_sigma_u":
            value = base * mp.e ** (-mp.mpf("0.5") * sigma * mp.mpf(args.omega) * u)
        elif variant == "apply_sigma_u":
            value = base * mp.e ** (mp.mpf("0.5") * sigma * mp.mpf(args.omega) * u)
        else:
            raise ValueError(variant)
        out[i] = value
    return out


def solve_feature_inverse(A: mp.matrix, y: mp.matrix, reg: mp.mpf):
    lhs = ctranspose(A) * A
    rhs = ctranspose(A) * y
    for i in range(min(lhs.rows, lhs.cols)):
        lhs[i, i] += reg
    return mp.lu_solve(lhs, rhs)


def solve_all_coefficients(args, xi: XiTransform, model: dict, variant: str):
    nodes = default_z_nodes()
    X = mp.matrix(args.basis, len(nodes))
    fit_errors = []
    for j, z in enumerate(nodes):
        y = target_vector(args, xi, z, model, variant)
        coeff = solve_feature_inverse(model["plusRows"], y, mp.mpf(args.reg))
        for i in range(args.basis):
            X[i, j] = coeff[i]
        fit_errors.append(to_float(frob_norm(model["plusRows"] * coeff - y) / max(frob_norm(y), mp.mpf("1e-300"))))
    return X, fit_errors


def joint_from_coefficients(model: dict, X: mp.matrix):
    P = herm(ctranspose(X) * model["P"] * X)
    M = herm(ctranspose(X) * model["M"] * X)
    C = ctranspose(X) * model["C"] * X
    return block2(P, C, ctranspose(C), M)


def fit_global_scale(mapped: mp.matrix, target: mp.matrix):
    num = mp.mpf("0")
    den = mp.mpf("0")
    for i in range(mapped.rows):
        for j in range(mapped.cols):
            num += mp.re(mp.conj(mapped[i, j]) * target[i, j])
            den += abs(mapped[i, j]) ** 2
    scale = max(mp.mpf("0"), num / den) if den else mp.mpf("0")
    return scale, rel_frob(scale * mapped, target)


def residual_report(mapped: mp.matrix, target: dict):
    n = target["P"].rows
    raw = {k: to_float(v) for k, v in block_residuals(mapped, target["joint"], n).items()}
    scale, scaled = fit_global_scale(mapped, target["joint"])
    raw["globalScale"] = to_float(scale)
    raw["scaledJointRelative"] = to_float(scaled)
    raw["mappedJointNorm"] = to_float(frob_norm(mapped))
    return raw


def trace_lift_coordinates(full_model: dict, trace_model: dict, coeffs: mp.matrix):
    X = mp.matrix(trace_model["traceRank"], coeffs.cols)
    for j in range(coeffs.cols):
        y = full_model["R"] * coeffs[:, j]
        c = solve_trace_coordinates(trace_model, y)
        for i in range(c.rows):
            X[i, j] = c[i]
    return X


def evaluate_variant(args, xi, full_model, trace_model, target, variant: str):
    coeffs, fit_errors = solve_all_coefficients(args, xi, full_model, variant)
    direct = joint_from_coefficients(full_model, coeffs)
    trace_coords = trace_lift_coordinates(full_model, trace_model, coeffs)
    lifted = mapped_joint(trace_model, trace_coords)
    return {
        "variant": variant,
        "featureFitRelativeMax": max(fit_errors),
        "featureFitRelativeMean": sum(fit_errors) / len(fit_errors),
        "coefficientNormFrobenius": to_float(frob_norm(coeffs)),
        "traceCoordinateNormFrobenius": to_float(frob_norm(trace_coords)),
        "directFeatureGram": residual_report(direct, target),
        "greenTraceGram": residual_report(lifted, target),
    }


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
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--json-out", default="klm_debranges_feature_inverse_candidate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    full_model = full_feature_model(args)
    trace_model = finite_trace_model(args)
    target = hardy_targets(xi, float(args.omega), None)
    variants = ["same", "split", "branch_unweighted", "undo_sigma_u", "apply_sigma_u"]
    results = [evaluate_variant(args, xi, full_model, trace_model, target, variant) for variant in variants]
    best_direct = min(results, key=lambda row: row["directFeatureGram"]["scaledJointRelative"])
    best_lifted = min(results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"])
    data = {
        "theoremName": "feature inverse candidate for primitive de Branges evaluation vector",
        "basis": args.basis,
        "constraints": args.constraints,
        "omega": float(args.omega),
        "traceRank": trace_model["traceRank"],
        "variants": variants,
        "results": results,
        "bestDirectFeatureGram": best_direct,
        "bestGreenTraceGram": best_lifted,
        "directFeatureMapConstructed": bool(best_direct["directFeatureGram"]["scaledJointRelative"] < 1e-6),
        "traceMapConstructed": bool(best_lifted["greenTraceGram"]["scaledJointRelative"] < 1e-6),
        "diagnosis": (
            "Finite inversion of the coefficient-level plus feature map was tested. "
            "If residuals remain large, the missing transform is not the naive "
            "row-grid inverse of G_+; it requires the exact continuous Weyl/"
            "Volterra inverse and branch normalization."
        ),
        "nextProofTarget": (
            "Derive the continuous adjoint normal equation G_+^*G_+ f_z=G_+^*h_z^+ "
            "and solve it analytically, including branch normalization and trace "
            "regularity, rather than relying on finite row-grid inversion."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Feature inverse candidate for primitive de Branges evaluation vector")
    print(f"  basis={args.basis} constraints={args.constraints} trace_rank={trace_model['traceRank']}")
    for row in results:
        print(
            f"  {row['variant']:17s} fit_max={row['featureFitRelativeMax']:.3e} "
            f"direct_scaled={row['directFeatureGram']['scaledJointRelative']:.3e} "
            f"trace_scaled={row['greenTraceGram']['scaledJointRelative']:.3e}"
        )
    print(
        f"  best direct={best_direct['variant']} "
        f"{best_direct['directFeatureGram']['scaledJointRelative']:.6e}"
    )
    print(
        f"  best trace={best_lifted['variant']} "
        f"{best_lifted['greenTraceGram']['scaledJointRelative']:.6e}"
    )
    print(f"  direct feature map constructed: {data['directFeatureMapConstructed']}")
    print(f"  trace map constructed: {data['traceMapConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
