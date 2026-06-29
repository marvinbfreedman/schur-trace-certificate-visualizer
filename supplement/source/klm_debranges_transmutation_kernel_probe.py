#!/usr/bin/env python3
r"""Probe Hardy-to-Volterra transmutation kernels.

The scalar and coupled normal equations failed when the Hardy feature was
identified with the Volterra output variable ``u``.  This script tests the next
structural possibility:

    r = s + u.

We lift Hardy half-line data h(r) to the Volterra triangle (s,u,sigma) through
simple kernels

    (U_theta h)(s,u,sigma) = theta(s,u,sigma) h(s+u),

then solve the coupled lifted least-squares problem for f_z.  The resulting
coefficients are still tested against the true integrated Volterra branch
Grams and against the Lambda trace map.

This is a search/proof-or-disproof artifact: if a natural theta makes the
joint Gram residual collapse, it identifies the missing transmutation.  If all
natural theta fail, the needed U is not a simple coarea multiplier on r=s+u.
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
    joint_from_forms,
    raw_lambda_traces,
    residual_report,
    trace_lift,
)
from klm_debranges_coupled_branch_normal_equation import PAIR_TRANSFORMS  # noqa: E402
from klm_debranges_pullback_probe import XiTransform  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    mapped_joint,
    to_float,
)
from quotient_factorization_mp import (  # noqa: E402
    endpoint_ratios,
    legendre_quadrature,
    poly_value,
    shifted_legendre_polys,
)
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


def a_value(s: mp.mpf, u: mp.mpf, pcs) -> mp.mpf:
    es = mp.e**s
    eu = mp.e**u
    total = mp.mpf("0")
    for ratio, beta, c in endpoint_ratios(s, pcs):
        total += ratio * mp.e ** (beta * u - c * es * (eu - 1))
    return total


def theta_value(kind: str, s: mp.mpf, u: mp.mpf, sigma: int, omega: mp.mpf, A: mp.mpf, eps: int, L: mp.mpf):
    r = s + u
    branch = mp.e ** (mp.mpf("0.5") * sigma * omega * r)
    if kind == "plain":
        return mp.mpf("1")
    if kind == "coarea":
        length = min(L, max(r, mp.mpf("1e-30")))
        return 1 / mp.sqrt(length)
    if kind == "branch_undo":
        return 1 / branch
    if kind == "A":
        return A
    if kind == "sqrtA":
        return mp.sqrt(abs(A))
    if kind == "A_branch":
        return A * branch
    if kind == "feature_weight":
        return A * branch * (1 + eps * r)
    if kind == "r_weight":
        return 1 + eps * r
    raise ValueError(kind)


def lifted_feature_rows(args: argparse.Namespace):
    omega = mp.mpf(args.omega)
    L = mp.mpf(args.L)
    pcs = pieces_for(args.kind)
    polys = shifted_legendre_polys(args.basis, L)
    s_pts, s_wts = legendre_quadrature(L, args.s_order)
    u_pts, u_wts = mp.gauss_quadrature(args.u_order, "legendre")
    # Finite u-window for the lifted probe; this is intentionally a search
    # model, while integrated Volterra Grams are still tested separately.
    umax = mp.mpf(args.u_max)
    u_pts = [mp.mpf("0.5") * umax * (x + 1) for x in u_pts]
    u_wts = [mp.mpf("0.5") * umax * w for w in u_wts]

    rows_plus = []
    rows_minus = []
    meta = []
    for sigma in branch_sigmas(omega):
        bw = branch_weight(omega)
        for iu, u in enumerate(u_pts):
            for is_, s in enumerate(s_pts):
                r = s + u
                A = a_value(s, u, pcs)
                root = mp.sqrt(s_wts[is_] * u_wts[iu] * bw / 4)
                branch = mp.e ** (mp.mpf("0.5") * sigma * omega * r)
                prow = []
                mrow = []
                for poly in polys:
                    base = root * A * branch * poly_value(poly, s)
                    prow.append(base * (1 + r))
                    mrow.append(base * (1 - r))
                rows_plus.append(prow)
                rows_minus.append(mrow)
                meta.append((s, u, sigma, A, root))
    return {
        "plus": mp.matrix(rows_plus),
        "minus": mp.matrix(rows_minus),
        "meta": meta,
        "polys": polys,
    }


def target_vector(args, xi, z: complex, lifted: dict, hardy_branch: str, phase, theta_kind: str, eps: int):
    omega = mp.mpf(args.omega)
    L = mp.mpf(args.L)
    amp = hardy_feature_scalar(xi, float(args.omega), z, hardy_branch)
    y = mp.matrix(len(lifted["meta"]), 1)
    for i, (s, u, sigma, A, root) in enumerate(lifted["meta"]):
        r = s + u
        theta = theta_value(theta_kind, s, u, sigma, omega, A, eps, L)
        # root without branch-feature factor makes target live in the same
        # lifted L2 quadrature as the feature rows.
        y[i] = root * mp.mpc(phase) * amp * mp.e ** (mp.j * mp.mpc(z) * r) * theta
    return y


def solve_for_transform(args, xi, lifted, pair_transform: str, theta_kind: str):
    nodes = default_z_nodes()
    A = ctranspose(lifted["plus"]) * lifted["plus"] + ctranspose(lifted["minus"]) * lifted["minus"]
    for i in range(A.rows):
        A[i, i] += mp.mpf(args.reg)
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    for j, z in enumerate(nodes):
        yp = target_vector(args, xi, z, lifted, plus_target[0], plus_target[1], theta_kind, 1)
        ym = target_vector(args, xi, z, lifted, minus_target[0], minus_target[1], theta_kind, -1)
        b = ctranspose(lifted["plus"]) * yp + ctranspose(lifted["minus"]) * ym
        coeff = mp.lu_solve(A, b)
        for k in range(args.basis):
            X[k, j] = coeff[k]
        res = A * coeff - b
        residuals.append(to_float(frob_norm(res) / max(frob_norm(b), mp.mpf("1e-300"))))
    return X, residuals


def lifted_joint(lifted: dict, X: mp.matrix):
    from klm_debranges_trace_map_constructor import block2, herm

    P = herm(ctranspose(X) * ctranspose(lifted["plus"]) * lifted["plus"] * X)
    M = herm(ctranspose(X) * ctranspose(lifted["minus"]) * lifted["minus"] * X)
    C = ctranspose(X) * ctranspose(lifted["plus"]) * lifted["minus"] * X
    return block2(P, C, ctranspose(C), M)


def evaluate(args, xi, forms, lifted, trace_model, target, pair_transform: str, theta_kind: str):
    X, residuals = solve_for_transform(args, xi, lifted, pair_transform, theta_kind)
    lifted_gram = lifted_joint(lifted, X)
    row = {
        "pairTransform": pair_transform,
        "theta": theta_kind,
        "normalResidualMax": max(residuals),
        "normalResidualMean": sum(residuals) / len(residuals),
        "coefficientNormFrobenius": to_float(frob_norm(X)),
        "liftedTriangleGram": residual_report(lifted_gram, target),
    }
    if forms is not None:
        direct = joint_from_forms(forms, X)
        row["integratedVolterraGram"] = residual_report(direct, target)
    else:
        row["integratedVolterraGram"] = None
    if trace_model is None and forms is not None:
        raw = raw_lambda_traces(forms, args, X)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["greenTraceGram"] = None
    elif trace_model is not None and forms is not None:
        raw, coords = trace_lift(forms, trace_model, args, X)
        lifted_trace = mapped_joint(trace_model, coords)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["traceCoordinateNormFrobenius"] = to_float(frob_norm(coords))
        row["greenTraceGram"] = residual_report(lifted_trace, target)
    else:
        row["greenTraceGram"] = None
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
    parser.add_argument("--s-order", type=int, default=10)
    parser.add_argument("--u-order", type=int, default=24)
    parser.add_argument("--u-max", default="10")
    parser.add_argument("--laguerre-order", type=int, default=32)
    parser.add_argument("--rhs-laguerre-order", type=int, default=40)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=700)
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--skip-trace-lift", action="store_true")
    parser.add_argument("--skip-integrated-forms", action="store_true")
    parser.add_argument("--thetas", default="plain,coarea,branch_undo,A,sqrtA,A_branch,r_weight")
    parser.add_argument("--pair-transforms", default="direct,swap")
    parser.add_argument("--json-out", default="klm_debranges_transmutation_kernel_probe.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    forms = None if args.skip_integrated_forms else assemble_feature_forms(args)
    lifted = lifted_feature_rows(args)
    trace_model = None
    if not args.skip_trace_lift and forms is not None:
        trace_model = finite_trace_model(
            SimpleNamespace(
                **{
                    **vars(args),
                    "u_order": max(40, args.rhs_laguerre_order),
                }
            )
        )
    target = hardy_targets(xi, float(args.omega), None)
    thetas = [item.strip() for item in args.thetas.split(",") if item.strip()]
    pair_transforms = [item.strip() for item in args.pair_transforms.split(",") if item.strip()]
    results = []
    for pair in pair_transforms:
        if pair not in PAIR_TRANSFORMS:
            raise ValueError(pair)
        for theta in thetas:
            results.append(evaluate(args, xi, forms, lifted, trace_model, target, pair, theta))
    integrated_results = [row for row in results if row["integratedVolterraGram"] is not None]
    best_integrated = (
        min(integrated_results, key=lambda row: row["integratedVolterraGram"]["scaledJointRelative"])
        if integrated_results
        else None
    )
    best_lifted = min(results, key=lambda row: row["liftedTriangleGram"]["scaledJointRelative"])
    trace_results = [row for row in results if row["greenTraceGram"] is not None]
    best_trace = min(trace_results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"]) if trace_results else None
    data = {
        "theoremName": "Hardy-to-Volterra transmutation kernel probe",
        "ansatz": "(U_theta h)(s,u,sigma)=theta(s,u,sigma) h(s+u)",
        "basis": args.basis,
        "constraints": args.constraints,
        "omega": float(args.omega),
        "sOrder": args.s_order,
        "uOrder": args.u_order,
        "uMax": float(args.u_max),
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "integratedFormsSkipped": bool(args.skip_integrated_forms),
        "thetas": thetas,
        "pairTransforms": pair_transforms,
        "results": results,
        "bestIntegratedVolterraGram": best_integrated,
        "bestLiftedTriangleGram": best_lifted,
        "bestGreenTraceGram": best_trace,
        "simpleCoareaTransmutationConstructed": bool(
            best_trace
            and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6
        ),
        "diagnosis": (
            "This tests the natural coarea ansatz r=s+u.  If residuals do not "
            "collapse, the missing U is not a scalar multiplier on the coarea "
            "surface; it must be a nonlocal integral transform in r or an "
            "operator-valued Volterra inverse."
        ),
        "nextProofTarget": (
            "Derive a nonlocal transmutation kernel U(r;s,u,sigma), likely from "
            "the Mellin/Laplace representation of Xi and the Volterra ratio "
            "A_s(u)=Psi(s+u)/Psi(s), rather than from scalar coarea weights."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Hardy-to-Volterra transmutation kernel probe")
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
            f"  {row['pairTransform']:6s} {row['theta']:13s} "
            f"integrated={row['integratedVolterraGram']['scaledJointRelative']:.3e}"
            if row["integratedVolterraGram"] is not None
            else f"  {row['pairTransform']:6s} {row['theta']:13s} integrated=skipped",
            end=" "
        )
        print(
            f"lifted={row['liftedTriangleGram']['scaledJointRelative']:.3e} "
            f"trace={trace_text}"
        )
    if best_integrated:
        print(
            f"  best integrated={best_integrated['pairTransform']}/{best_integrated['theta']} "
            f"{best_integrated['integratedVolterraGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best integrated=skipped")
    print(
        f"  best lifted={best_lifted['pairTransform']}/{best_lifted['theta']} "
        f"{best_lifted['liftedTriangleGram']['scaledJointRelative']:.6e}"
    )
    if best_trace:
        print(
            f"  best trace={best_trace['pairTransform']}/{best_trace['theta']} "
            f"{best_trace['greenTraceGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best trace=skipped")
    print(f"  simple coarea transmutation constructed: {data['simpleCoareaTransmutationConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
