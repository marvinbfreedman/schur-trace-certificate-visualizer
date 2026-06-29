#!/usr/bin/env python3
r"""Coupled signed-branch normal equation for the KLM/de Branges bridge.

The plus-only equation

    G_+^*G_+ f_z = G_+^*h_z^+

does not recover the joint Hardy branch Grams.  This script tests the next
natural system: fit both signed branches at once,

    (G_+^*G_+ + G_-^*G_-) f_z
      = G_+^*T_+h_z + G_-^*T_-h_z,

where the finite list of ``pair transforms`` tries the obvious missing
Hardy-to-Volterra normalizations: direct branch matching, swapped branches,
minus sign, and phase rotations of the minus branch.  The feature kernels are
continuous Weyl/Volterra kernels; quadrature is used only to assemble their
Galerkin entries.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import default_z_nodes, hardy_feature_scalar  # noqa: E402
from klm_debranges_continuous_normal_equation import (  # noqa: E402
    assemble_feature_forms,
    branch_sigmas,
    branch_weight,
    fit_global_scale,
    one_feature_integral,
    raw_lambda_traces,
    residual_report,
    solve_normal,
    target_multiplier,
    trace_lift,
)
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    herm,
    mapped_joint,
    to_float,
)
from quotient_factorization_mp import endpoint_ratios, poly_value  # noqa: E402
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


PAIR_TRANSFORMS = {
    "direct": (("plus", 1), ("minus", 1)),
    "minus_neg": (("plus", 1), ("minus", -1)),
    "minus_i": (("plus", 1), ("minus", 1j)),
    "minus_minus_i": (("plus", 1), ("minus", -1j)),
    "swap": (("minus", 1), ("plus", 1)),
    "swap_minus_neg": (("minus", 1), ("plus", -1)),
}


def rhs_density_at_s_general(
    args: argparse.Namespace,
    xi: XiTransform,
    z: complex,
    s: mp.mpf,
    pcs,
    lag_nodes,
    lag_weights,
    variant: str,
    feature_eps: int,
    hardy_branch: str,
    phase,
):
    omega = mp.mpf(args.omega)
    amp = hardy_feature_scalar(xi, float(args.omega), z, hardy_branch)
    total = mp.mpc("0")
    w = branch_weight(omega)
    es = mp.e**s
    for sigma in branch_sigmas(omega):
        const, adjust = target_multiplier(args, sigma, mp.j * mp.mpc(z), variant)
        sigma_factor = mp.e ** (mp.mpf("0.5") * sigma * omega * s)
        for ratio_i, beta_i, c_i in endpoint_ratios(s, pcs):
            alpha = c_i * es
            p = beta_i + mp.mpf("0.5") * sigma * omega + adjust - 1
            total += (
                mp.sqrt(w / 4)
                * const
                * mp.mpc(phase)
                * amp
                * ratio_i
                * sigma_factor
                * one_feature_integral(alpha, p, s, feature_eps, lag_nodes, lag_weights)
            )
    return total


def assemble_branch_rhs(
    args: argparse.Namespace,
    xi: XiTransform,
    forms: dict,
    z: complex,
    variant: str,
    feature_eps: int,
    hardy_branch: str,
    phase,
):
    pcs = pieces_for(args.kind)
    lag_nodes, lag_weights = mp.gauss_quadrature(args.rhs_laguerre_order, "laguerre")
    b = mp.matrix(args.basis, 1)
    for ia, s in enumerate(forms["sPts"]):
        density = rhs_density_at_s_general(
            args, xi, z, s, pcs, lag_nodes, lag_weights, variant, feature_eps, hardy_branch, phase
        )
        for i, poly in enumerate(forms["polys"]):
            b[i] += forms["sWts"][ia] * poly_value(poly, s) * density
    return b


def pair_rhs(args, xi, forms, z: complex, variant: str, pair_transform: str):
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    b_plus = assemble_branch_rhs(args, xi, forms, z, variant, 1, plus_target[0], plus_target[1])
    b_minus = assemble_branch_rhs(args, xi, forms, z, variant, -1, minus_target[0], minus_target[1])
    return b_plus + b_minus


def build_rhs_cache(args, xi, forms, variants: list[str]):
    nodes = default_z_nodes()
    cache = {}
    for variant in variants:
        for j, z in enumerate(nodes):
            for feature_eps in (1, -1):
                for hardy_branch in ("plus", "minus"):
                    cache[(variant, j, feature_eps, hardy_branch)] = assemble_branch_rhs(
                        args, xi, forms, z, variant, feature_eps, hardy_branch, 1
                    )
    return cache


def pair_rhs_cached(cache: dict, variant: str, node_index: int, pair_transform: str):
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    return (
        mp.mpc(plus_target[1]) * cache[(variant, node_index, 1, plus_target[0])]
        + mp.mpc(minus_target[1]) * cache[(variant, node_index, -1, minus_target[0])]
    )


def joint_from_forms(forms: dict, X: mp.matrix):
    P = herm(ctranspose(X) * forms["P"] * X)
    M = herm(ctranspose(X) * forms["M"] * X)
    C = ctranspose(X) * forms["C"] * X
    from klm_debranges_trace_map_constructor import block2  # local import avoids export assumptions

    return block2(P, C, ctranspose(C), M)


def solve_coupled_coefficients(args, xi, forms, variant: str, pair_transform: str, rhs_cache: dict | None):
    nodes = default_z_nodes()
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    A = herm(forms["P"] + forms["M"])
    for j, z in enumerate(nodes):
        b = (
            pair_rhs_cached(rhs_cache, variant, j, pair_transform)
            if rhs_cache is not None
            else pair_rhs(args, xi, forms, z, variant, pair_transform)
        )
        coeff = solve_normal(A, b, mp.mpf(args.reg))
        for i in range(args.basis):
            X[i, j] = coeff[i]
        res = A * coeff - b
        residuals.append(to_float(frob_norm(res) / max(frob_norm(b), mp.mpf("1e-300"))))
    return X, residuals


def evaluate(args, xi, forms, trace_model, target, variant: str, pair_transform: str, rhs_cache: dict | None):
    X, residuals = solve_coupled_coefficients(args, xi, forms, variant, pair_transform, rhs_cache)
    direct = joint_from_forms(forms, X)
    row = {
        "variant": variant,
        "pairTransform": pair_transform,
        "normalEquationResidualMax": max(residuals),
        "normalEquationResidualMean": sum(residuals) / len(residuals),
        "coefficientNormFrobenius": to_float(frob_norm(X)),
        "directContinuousGram": residual_report(direct, target),
    }
    if trace_model is None:
        raw = raw_lambda_traces(forms, args, X)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["greenTraceGram"] = None
    else:
        raw, coords = trace_lift(forms, trace_model, args, X)
        lifted = mapped_joint(trace_model, coords)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["traceCoordinateNormFrobenius"] = to_float(frob_norm(coords))
        row["greenTraceGram"] = residual_report(lifted, target)
    return row


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
    parser.add_argument("--s-order", type=int, default=16)
    parser.add_argument("--laguerre-order", type=int, default=40)
    parser.add_argument("--rhs-laguerre-order", type=int, default=60)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=800)
    parser.add_argument("--hardy-rmax", type=float, default=40.0)
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--skip-trace-lift", action="store_true")
    parser.add_argument("--variants", default="same,undo_sigma_u,apply_sigma_u")
    parser.add_argument("--pair-transforms", default=",".join(PAIR_TRANSFORMS))
    parser.add_argument("--json-out", default="klm_debranges_coupled_branch_normal_equation.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    forms = assemble_feature_forms(args)
    trace_model = None
    if not args.skip_trace_lift:
        trace_model = finite_trace_model(
            SimpleNamespace(
                **{
                    **vars(args),
                    "u_order": max(40, args.rhs_laguerre_order),
                    "u_max": "10",
                }
            )
        )
    target = hardy_targets(xi, float(args.omega), None)
    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    pair_transforms = [item.strip() for item in args.pair_transforms.split(",") if item.strip()]
    for item in pair_transforms:
        if item not in PAIR_TRANSFORMS:
            raise ValueError(f"unknown pair transform {item}")
    results = []
    rhs_cache = build_rhs_cache(args, xi, forms, variants)
    for pair_transform in pair_transforms:
        for variant in variants:
            results.append(evaluate(args, xi, forms, trace_model, target, variant, pair_transform, rhs_cache))
    best_direct = min(results, key=lambda row: row["directContinuousGram"]["scaledJointRelative"])
    trace_results = [row for row in results if row["greenTraceGram"] is not None]
    best_trace = min(trace_results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"]) if trace_results else None
    data = {
        "theoremName": "coupled signed-branch Weyl/Volterra normal equation",
        "normalEquation": "(G_+^*G_+ + G_-^*G_-)f_z = G_+^*T_+h_z + G_-^*T_-h_z",
        "basis": args.basis,
        "constraints": args.constraints,
        "omega": float(args.omega),
        "sOrder": args.s_order,
        "laguerreOrder": args.laguerre_order,
        "rhsLaguerreOrder": args.rhs_laguerre_order,
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "variants": variants,
        "pairTransforms": pair_transforms,
        "results": results,
        "bestDirectContinuousGram": best_direct,
        "bestGreenTraceGram": best_trace,
        "coupledNormalEquationSolved": True,
        "directContinuousMapConstructed": bool(best_direct["directContinuousGram"]["scaledJointRelative"] < 1e-6),
        "traceMapConstructed": bool(best_trace and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6),
        "diagnosis": (
            "The coupled signed-branch normal equations test the obvious missing "
            "normalizations: direct, swapped, and phase-rotated branch targets. "
            "If residuals remain large, the missing transform is not a scalar "
            "branch renormalization; it must mix the Hardy half-line variable with "
            "the Volterra endpoint variable more substantially."
        ),
        "nextProofTarget": (
            "If the coupled branch system fails, derive an integral "
            "Hardy-to-Volterra transmutation kernel U(r;s,u,sigma) rather than "
            "a scalar branch normalization."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Coupled signed-branch Weyl/Volterra normal equation")
    print(
        f"  basis={args.basis} constraints={args.constraints} "
        f"trace_rank={trace_model['traceRank'] if trace_model is not None else 'skipped'}"
    )
    for row in results:
        trace_text = (
            f"{row['greenTraceGram']['scaledJointRelative']:.3e}"
            if row["greenTraceGram"] is not None
            else "skipped"
        )
        print(
            f"  {row['pairTransform']:14s} {row['variant']:13s} "
            f"res={row['normalEquationResidualMax']:.1e} "
            f"direct={row['directContinuousGram']['scaledJointRelative']:.3e} "
            f"trace={trace_text}"
        )
    print(
        f"  best direct={best_direct['pairTransform']}/{best_direct['variant']} "
        f"{best_direct['directContinuousGram']['scaledJointRelative']:.6e}"
    )
    if best_trace:
        print(
            f"  best trace={best_trace['pairTransform']}/{best_trace['variant']} "
            f"{best_trace['greenTraceGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best trace=skipped")
    print(f"  direct continuous map constructed: {data['directContinuousMapConstructed']}")
    print(f"  trace map constructed: {data['traceMapConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
