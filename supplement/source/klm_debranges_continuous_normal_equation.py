#!/usr/bin/env python3
r"""Continuous Galerkin normal equation for the de Branges trace map.

This is the next version after ``klm_debranges_feature_inverse_candidate.py``.
Instead of solving a sampled row equation

    G_+ f_z ~= h_z^+

we assemble the continuous Galerkin normal equation

    G_+^* G_+ f_z = G_+^* h_z^+

from the Weyl/Volterra feature kernel itself.  The only numerical
discretization is quadrature of the displayed continuous integral formulas.
After solving, we compute

    x_z(a_j) = Lambda_{a_j}(f_z).

For a sign branch sigma and r=s+u, the feature factors are

    g_+(s,u,sigma) = sqrt(w_sigma/4) A_s(u) e^{sigma omega r/2} (1+r),
    g_-(s,u,sigma) = sqrt(w_sigma/4) A_s(u) e^{sigma omega r/2} (1-r),

where A_s(u)=Psi(s+u)/Psi(s).  This script builds the continuous plus,
minus, and cross Gram matrices induced by those features.
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
from klm_debranges_feature_inverse_candidate import fit_global_scale, residual_report  # noqa: E402
from klm_debranges_lambda_trace_candidate import solve_trace_coordinates  # noqa: E402
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    block2,
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    herm,
    local_args,
    mapped_joint,
    sym,
    to_float,
)
from quotient_factorization_mp import (  # noqa: E402
    endpoint_ratios,
    legendre_quadrature,
    poly_value,
    shifted_legendre_polys,
    trace_matrix,
)
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


def branch_sigmas(omega: mp.mpf) -> list[int]:
    return [0] if omega == 0 else [-1, 1]


def branch_weight(omega: mp.mpf) -> mp.mpf:
    return mp.mpf("1") if omega == 0 else mp.mpf("0.5")


def q_power_integral(alpha: mp.mpf, p, center_s: mp.mpf, center_t: mp.mpf, eps_l: int, eps_r: int, nodes, weights):
    """Integral of (1+eps_l(s+u))(1+eps_r(t+u)) q^p e^{-alpha(q-1)} du."""
    total = mp.mpc("0")
    for node, weight in zip(nodes, weights):
        q = 1 + node / alpha
        u = mp.log(q)
        left = 1 + eps_l * (center_s + u)
        right = 1 + eps_r * (center_t + u)
        total += weight * left * right * (q ** p)
    return total / alpha


def one_feature_integral(alpha: mp.mpf, p, s: mp.mpf, eps: int, nodes, weights):
    """Integral of (1+eps(s+u)) q^p e^{-alpha(q-1)} du."""
    total = mp.mpc("0")
    for node, weight in zip(nodes, weights):
        q = 1 + node / alpha
        u = mp.log(q)
        total += weight * (1 + eps * (s + u)) * (q ** p)
    return total / alpha


def feature_kernel(s: mp.mpf, t: mp.mpf, omega: mp.mpf, pcs, eps_l: int, eps_r: int, lag_nodes, lag_weights):
    left = endpoint_ratios(s, pcs)
    right = endpoint_ratios(t, pcs)
    es = mp.e**s
    et = mp.e**t
    total = mp.mpc("0")
    w = branch_weight(omega)
    for sigma in branch_sigmas(omega):
        sigma_factor = mp.e ** (mp.mpf("0.5") * sigma * omega * (s + t))
        for ratio_i, beta_i, c_i in left:
            for ratio_j, beta_j, c_j in right:
                alpha = c_i * es + c_j * et
                p = beta_i + beta_j + sigma * omega - 1
                total += (
                    (w / 4)
                    * ratio_i
                    * ratio_j
                    * sigma_factor
                    * q_power_integral(alpha, p, s, t, eps_l, eps_r, lag_nodes, lag_weights)
                )
    return total


def assemble_feature_forms(args: argparse.Namespace):
    omega = mp.mpf(args.omega)
    L = mp.mpf(args.L)
    pcs = pieces_for(args.kind)
    polys = shifted_legendre_polys(args.basis, L)
    s_pts, s_wts = legendre_quadrature(L, args.s_order)
    lag_nodes, lag_weights = mp.gauss_quadrature(args.laguerre_order, "laguerre")
    values = mp.matrix(len(s_pts), args.basis)
    for i, s in enumerate(s_pts):
        for j, poly in enumerate(polys):
            values[i, j] = poly_value(poly, s)

    P = mp.matrix(args.basis)
    M = mp.matrix(args.basis)
    C = mp.matrix(args.basis)
    for ia, s in enumerate(s_pts):
        for ib, t in enumerate(s_pts):
            weight = s_wts[ia] * s_wts[ib]
            kpp = feature_kernel(s, t, omega, pcs, 1, 1, lag_nodes, lag_weights)
            kmm = feature_kernel(s, t, omega, pcs, -1, -1, lag_nodes, lag_weights)
            kpm = feature_kernel(s, t, omega, pcs, 1, -1, lag_nodes, lag_weights)
            for i in range(args.basis):
                vi = values[ia, i]
                for j in range(args.basis):
                    vj = values[ib, j]
                    fac = weight * vi * vj
                    P[i, j] += fac * kpp
                    M[i, j] += fac * kmm
                    C[i, j] += fac * kpm
    return {
        "P": herm(P),
        "M": herm(M),
        "C": C,
        "polys": polys,
        "sPts": s_pts,
        "sWts": s_wts,
        "laguerreOrder": args.laguerre_order,
        "sOrder": args.s_order,
    }


def target_multiplier(args: argparse.Namespace, sigma: int, u_power_adjust: mp.mpf, variant: str):
    if variant == "same":
        return mp.mpf("1"), u_power_adjust
    if variant == "split":
        return (mp.mpf("1") / mp.sqrt(2) if sigma != 0 else mp.mpf("1")), u_power_adjust
    if variant == "branch_unweighted":
        return mp.mpf("1") / mp.sqrt(branch_weight(mp.mpf(args.omega))), u_power_adjust
    if variant == "undo_sigma_u":
        return mp.mpf("1"), u_power_adjust - mp.mpf("0.5") * sigma * mp.mpf(args.omega)
    if variant == "apply_sigma_u":
        return mp.mpf("1"), u_power_adjust + mp.mpf("0.5") * sigma * mp.mpf(args.omega)
    raise ValueError(variant)


def rhs_density_at_s(args: argparse.Namespace, xi: XiTransform, z: complex, s: mp.mpf, pcs, lag_nodes, lag_weights, variant: str):
    omega = mp.mpf(args.omega)
    amp = hardy_feature_scalar(xi, float(args.omega), z, "plus")
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
                * amp
                * ratio_i
                * sigma_factor
                * one_feature_integral(alpha, p, s, 1, lag_nodes, lag_weights)
            )
    return total


def assemble_rhs(args: argparse.Namespace, xi: XiTransform, forms: dict, z: complex, variant: str):
    pcs = pieces_for(args.kind)
    lag_nodes, lag_weights = mp.gauss_quadrature(args.rhs_laguerre_order, "laguerre")
    b = mp.matrix(args.basis, 1)
    for ia, s in enumerate(forms["sPts"]):
        density = rhs_density_at_s(args, xi, z, s, pcs, lag_nodes, lag_weights, variant)
        for i, poly in enumerate(forms["polys"]):
            b[i] += forms["sWts"][ia] * poly_value(poly, s) * density
    return b


def solve_normal(P: mp.matrix, b: mp.matrix, reg: mp.mpf):
    lhs = mp.matrix(P)
    for i in range(lhs.rows):
        lhs[i, i] += reg
    return mp.lu_solve(lhs, b)


def solve_coefficients(args: argparse.Namespace, xi: XiTransform, forms: dict, variant: str):
    nodes = default_z_nodes()
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    for j, z in enumerate(nodes):
        b = assemble_rhs(args, xi, forms, z, variant)
        coeff = solve_normal(forms["P"], b, mp.mpf(args.reg))
        for i in range(args.basis):
            X[i, j] = coeff[i]
        res = forms["P"] * coeff - b
        residuals.append(to_float(frob_norm(res) / max(frob_norm(b), mp.mpf("1e-300"))))
    return X, residuals


def joint_from_forms(forms: dict, X: mp.matrix):
    P = herm(ctranspose(X) * forms["P"] * X)
    M = herm(ctranspose(X) * forms["M"] * X)
    C = ctranspose(X) * forms["C"] * X
    return block2(P, C, ctranspose(C), M)


def raw_lambda_traces(forms: dict, args: argparse.Namespace, X: mp.matrix):
    _centers, R = trace_matrix(forms["polys"], local_args(args))
    raw_traces = mp.matrix(R.rows, X.cols)
    for j in range(X.cols):
        y = R * X[:, j]
        for i in range(y.rows):
            raw_traces[i, j] = y[i]
    return raw_traces


def trace_lift(forms: dict, trace_model: dict, args: argparse.Namespace, X: mp.matrix):
    _centers, R = trace_matrix(forms["polys"], local_args(args))
    trace_coords = mp.matrix(trace_model["traceRank"], X.cols)
    raw_traces = mp.matrix(R.rows, X.cols)
    for j in range(X.cols):
        y = R * X[:, j]
        for i in range(y.rows):
            raw_traces[i, j] = y[i]
        c = solve_trace_coordinates(trace_model, y)
        for i in range(c.rows):
            trace_coords[i, j] = c[i]
    return raw_traces, trace_coords


def evaluate_variant(args, xi, forms, trace_model, target, variant: str):
    X, normal_residuals = solve_coefficients(args, xi, forms, variant)
    direct = joint_from_forms(forms, X)
    row = {
        "variant": variant,
        "normalEquationResidualMax": max(normal_residuals),
        "normalEquationResidualMean": sum(normal_residuals) / len(normal_residuals),
        "coefficientNormFrobenius": to_float(frob_norm(X)),
        "directContinuousGram": residual_report(direct, target),
    }
    if trace_model is None:
        raw_traces = raw_lambda_traces(forms, args, X)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw_traces))
        row["greenTraceGram"] = None
    else:
        raw_traces, trace_coords = trace_lift(forms, trace_model, args, X)
        lifted = mapped_joint(trace_model, trace_coords)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw_traces))
        row["traceCoordinateNormFrobenius"] = to_float(frob_norm(trace_coords))
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
    parser.add_argument("--s-order", type=int, default=18)
    parser.add_argument("--laguerre-order", type=int, default=50)
    parser.add_argument("--rhs-laguerre-order", type=int, default=70)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=800)
    parser.add_argument("--hardy-rmax", type=float, default=40.0)
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--skip-trace-lift", action="store_true")
    parser.add_argument("--json-out", default="klm_debranges_continuous_normal_equation.json")
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
    variants = ["same", "split", "branch_unweighted", "undo_sigma_u", "apply_sigma_u"]
    results = [evaluate_variant(args, xi, forms, trace_model, target, variant) for variant in variants]
    best_direct = min(results, key=lambda row: row["directContinuousGram"]["scaledJointRelative"])
    trace_results = [row for row in results if row["greenTraceGram"] is not None]
    best_trace = min(trace_results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"]) if trace_results else None
    data = {
        "theoremName": "continuous Weyl/Volterra normal equation for de Branges evaluation trace",
        "normalEquation": "G_+^*G_+ f_z = G_+^*h_z^+",
        "plusKernel": (
            "N_+(s,t)=sum_sigma w_sigma/4 int_0^inf "
            "(1+s+u)(1+t+u) A_s(u)A_t(u) "
            "exp(sigma omega (u+(s+t)/2)) du"
        ),
        "rhsKernel": (
            "b_z(s)=sum_sigma sqrt(w_sigma/4) int_0^inf "
            "(1+s+u) A_s(u) exp(sigma omega(s+u)/2) h_z^+(u,sigma) du"
        ),
        "basis": args.basis,
        "constraints": args.constraints,
        "omega": float(args.omega),
        "sOrder": args.s_order,
        "laguerreOrder": args.laguerre_order,
        "rhsLaguerreOrder": args.rhs_laguerre_order,
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "variants": variants,
        "results": results,
        "bestDirectContinuousGram": best_direct,
        "bestGreenTraceGram": best_trace,
        "continuousNormalEquationSolved": True,
        "directContinuousMapConstructed": bool(best_direct["directContinuousGram"]["scaledJointRelative"] < 1e-6),
        "traceMapConstructed": bool(best_trace and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6),
        "diagnosis": (
            "The continuous Galerkin normal equation has been assembled from the "
            "Weyl/Volterra kernel and solved.  If the residuals remain large under "
            "refinement, the missing bridge is not this plus-branch least-squares "
            "inverse; the target must include an additional transform or a different "
            "Hardy-to-Volterra branch normalization."
        ),
        "nextProofTarget": (
            "Refine the continuous normal-equation certificate and derive the "
            "missing Hardy-to-Volterra branch normalization or extra transform "
            "if the trace-lifted joint Gram remains far from the Hardy joint Gram."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Continuous Weyl/Volterra normal equation")
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
            f"  {row['variant']:17s} normal_res={row['normalEquationResidualMax']:.3e} "
            f"direct_scaled={row['directContinuousGram']['scaledJointRelative']:.3e} "
            f"trace_scaled={trace_text}"
        )
    print(
        f"  best direct={best_direct['variant']} "
        f"{best_direct['directContinuousGram']['scaledJointRelative']:.6e}"
    )
    if best_trace is not None:
        print(
            f"  best trace={best_trace['variant']} "
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
