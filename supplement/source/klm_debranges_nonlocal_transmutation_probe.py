#!/usr/bin/env python3
r"""Nonlocal Hardy-to-Volterra transmutation kernel probe.

The scalar coarea map

    (U_theta h)(s,u,sigma)=theta(s,u,sigma) h(s+u)

did not close the KLM/de Branges bridge.  This script tests genuinely
nonlocal half-line kernels.  For h_z(r)=E(z)e^{izr}, the tested kernels have
explicit Fourier actions, for example

    int_c^inf exp(-lambda(r-c)) e^{izr} dr = exp(izc)/(lambda-iz),

and the half-line two-sided resolvent centered at c adds the backward term

    int_0^c exp(-lambda(c-r)) e^{izr} dr
      = (exp(izc)-exp(-lambda c))/(lambda+iz).

The centers c are chosen from the Volterra/Mellin geometry:

    c=s+u,
    c=u,
    c=log(1+exp(s)(exp(u)-1)).

This is still a probe, not a theorem: it checks whether a small explicit
Mellin/Laplace dictionary contains the missing transmutation.
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
from klm_debranges_transmutation_kernel_probe import lifted_feature_rows, theta_value  # noqa: E402


def center_value(kind: str, s: mp.mpf, u: mp.mpf) -> mp.mpf:
    if kind == "su":
        return s + u
    if kind == "u":
        return u
    if kind == "mellin_eta":
        return mp.log(1 + mp.e**s * (mp.e**u - 1))
    if kind == "avg_su_eta":
        return mp.mpf("0.5") * ((s + u) + mp.log(1 + mp.e**s * (mp.e**u - 1)))
    raise ValueError(kind)


def fourier_kernel(kind: str, z: complex, center: mp.mpf, lam: mp.mpf):
    zz = mp.mpc(z)
    iz = mp.j * zz
    if kind == "coarea":
        return mp.e ** (iz * center)
    if kind == "forward":
        return mp.e ** (iz * center) / (lam - iz)
    if kind == "backward":
        return (mp.e ** (iz * center) - mp.e ** (-lam * center)) / (lam + iz)
    if kind == "twosided":
        return (
            mp.e ** (iz * center) / (lam - iz)
            + (mp.e ** (iz * center) - mp.e ** (-lam * center)) / (lam + iz)
        )
    if kind == "forward2":
        return mp.e ** (iz * center) / ((lam - iz) ** 2)
    if kind == "backward2":
        # int_0^c (c-r) exp(-lambda(c-r)) exp(izr) dr.
        a = lam + iz
        return (mp.e ** (iz * center) * (a * center - 1) + mp.e ** (-lam * center)) / (a**2)
    if kind == "heat":
        tau = 1 / max(lam, mp.mpf("1e-30"))
        return mp.e ** (iz * center - mp.mpf("0.5") * tau * zz * zz)
    raise ValueError(kind)


def target_vector(args, xi, z: complex, lifted: dict, hardy_branch: str, phase, spec: dict, eps: int):
    omega = mp.mpf(args.omega)
    L = mp.mpf(args.L)
    amp = hardy_feature_scalar(xi, float(args.omega), z, hardy_branch)
    y = mp.matrix(len(lifted["meta"]), 1)
    lam = mp.mpf(spec["lambda"])
    for i, (s, u, sigma, A, root) in enumerate(lifted["meta"]):
        c = center_value(spec["center"], s, u)
        theta = theta_value(spec["theta"], s, u, sigma, omega, A, eps, L)
        y[i] = root * mp.mpc(phase) * amp * theta * fourier_kernel(spec["kernel"], z, c, lam)
    return y


def solve_for_spec(args, xi, lifted, pair_transform: str, spec: dict):
    nodes = default_z_nodes()
    A = ctranspose(lifted["plus"]) * lifted["plus"] + ctranspose(lifted["minus"]) * lifted["minus"]
    for i in range(A.rows):
        A[i, i] += mp.mpf(args.reg)
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    for j, z in enumerate(nodes):
        yp = target_vector(args, xi, z, lifted, plus_target[0], plus_target[1], spec, 1)
        ym = target_vector(args, xi, z, lifted, minus_target[0], minus_target[1], spec, -1)
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


def evaluate(args, xi, forms, lifted, trace_model, target, pair_transform: str, spec: dict):
    X, residuals = solve_for_spec(args, xi, lifted, pair_transform, spec)
    lift = lifted_joint(lifted, X)
    row = {
        "pairTransform": pair_transform,
        "spec": spec,
        "normalResidualMax": max(residuals),
        "normalResidualMean": sum(residuals) / len(residuals),
        "coefficientNormFrobenius": to_float(frob_norm(X)),
        "liftedTriangleGram": residual_report(lift, target),
    }
    if forms is not None:
        direct = joint_from_forms(forms, X)
        row["integratedVolterraGram"] = residual_report(direct, target)
    else:
        row["integratedVolterraGram"] = None
    if forms is not None and trace_model is None:
        raw = raw_lambda_traces(forms, args, X)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["greenTraceGram"] = None
    elif forms is not None and trace_model is not None:
        raw, coords = trace_lift(forms, trace_model, args, X)
        lifted_trace = mapped_joint(trace_model, coords)
        row["rawTraceNormFrobenius"] = to_float(frob_norm(raw))
        row["traceCoordinateNormFrobenius"] = to_float(frob_norm(coords))
        row["greenTraceGram"] = residual_report(lifted_trace, target)
    else:
        row["greenTraceGram"] = None
    return row


def parse_specs(args: argparse.Namespace):
    kernels = [x.strip() for x in args.kernels.split(",") if x.strip()]
    centers = [x.strip() for x in args.centers.split(",") if x.strip()]
    thetas = [x.strip() for x in args.thetas.split(",") if x.strip()]
    lambdas = [x.strip() for x in args.lambdas.split(",") if x.strip()]
    specs = []
    for kernel in kernels:
        for center in centers:
            for theta in thetas:
                for lam in lambdas:
                    specs.append({"kernel": kernel, "center": center, "theta": theta, "lambda": lam})
    return specs


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
    parser.add_argument("--kernels", default="forward,backward,twosided,forward2,heat")
    parser.add_argument("--centers", default="su,u,mellin_eta")
    parser.add_argument("--thetas", default="plain,coarea,A,sqrtA")
    parser.add_argument("--lambdas", default="0.5,1,2,4")
    parser.add_argument("--pair-transforms", default="direct,swap")
    parser.add_argument("--max-specs", type=int, default=0)
    parser.add_argument("--json-out", default="klm_debranges_nonlocal_transmutation_probe.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    forms = None if args.skip_integrated_forms else assemble_feature_forms(args)
    lifted = lifted_feature_rows(args)
    trace_model = None
    if forms is not None and not args.skip_trace_lift:
        trace_model = finite_trace_model(
            SimpleNamespace(
                **{
                    **vars(args),
                    "u_order": max(40, args.rhs_laguerre_order),
                }
            )
        )
    target = hardy_targets(xi, float(args.omega), None)
    pair_transforms = [x.strip() for x in args.pair_transforms.split(",") if x.strip()]
    specs = parse_specs(args)
    if args.max_specs > 0:
        specs = specs[: args.max_specs]
    results = []
    for pair in pair_transforms:
        if pair not in PAIR_TRANSFORMS:
            raise ValueError(pair)
        for spec in specs:
            results.append(evaluate(args, xi, forms, lifted, trace_model, target, pair, spec))
    best_lifted = min(results, key=lambda row: row["liftedTriangleGram"]["scaledJointRelative"])
    integrated_results = [row for row in results if row["integratedVolterraGram"] is not None]
    best_integrated = (
        min(integrated_results, key=lambda row: row["integratedVolterraGram"]["scaledJointRelative"])
        if integrated_results
        else None
    )
    trace_results = [row for row in results if row["greenTraceGram"] is not None]
    best_trace = min(trace_results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"]) if trace_results else None
    data = {
        "theoremName": "nonlocal Hardy-to-Volterra transmutation kernel probe",
        "ansatz": "U generated by explicit half-line resolvent/heat kernels centered at Volterra/Mellin locations",
        "basis": args.basis,
        "constraints": args.constraints,
        "omega": float(args.omega),
        "sOrder": args.s_order,
        "uOrder": args.u_order,
        "uMax": float(args.u_max),
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "integratedFormsSkipped": bool(args.skip_integrated_forms),
        "pairTransforms": pair_transforms,
        "specCount": len(specs),
        "results": results,
        "bestLiftedTriangleGram": best_lifted,
        "bestIntegratedVolterraGram": best_integrated,
        "bestGreenTraceGram": best_trace,
        "nonlocalDictionaryConstructedBridge": bool(
            best_trace and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6
        ),
        "diagnosis": (
            "This tests a finite Mellin/Laplace dictionary of nonlocal kernels. "
            "If residuals remain large, the missing U is not a simple resolvent "
            "or heat kernel centered at s+u, u, or log(1+e^s(e^u-1))."
        ),
        "nextProofTarget": (
            "Derive U from the exact Xi Mellin transform by matching the Volterra "
            "mode exponentials term-by-term, rather than using a generic "
            "resolvent/heat dictionary."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Nonlocal Hardy-to-Volterra transmutation probe")
    print(
        f"  basis={args.basis} constraints={args.constraints} specs={len(specs)} "
        f"trace_rank={trace_model['traceRank'] if trace_model else 'skipped'}"
    )
    for row in results:
        spec = row["spec"]
        trace_text = (
            f"{row['greenTraceGram']['scaledJointRelative']:.3e}"
            if row["greenTraceGram"] is not None
            else "skipped"
        )
        integ_text = (
            f"{row['integratedVolterraGram']['scaledJointRelative']:.3e}"
            if row["integratedVolterraGram"] is not None
            else "skipped"
        )
        print(
            f"  {row['pairTransform']:6s} {spec['kernel']:9s} {spec['center']:10s} "
            f"{spec['theta']:8s} lam={spec['lambda']:>4s} "
            f"lifted={row['liftedTriangleGram']['scaledJointRelative']:.3e} "
            f"integrated={integ_text} trace={trace_text}"
        )
    b = best_lifted
    print(
        f"  best lifted={b['pairTransform']}/{b['spec']['kernel']}/"
        f"{b['spec']['center']}/{b['spec']['theta']}/lambda={b['spec']['lambda']} "
        f"{b['liftedTriangleGram']['scaledJointRelative']:.6e}"
    )
    if best_integrated:
        b = best_integrated
        print(
            f"  best integrated={b['pairTransform']}/{b['spec']['kernel']}/"
            f"{b['spec']['center']}/{b['spec']['theta']}/lambda={b['spec']['lambda']} "
            f"{b['integratedVolterraGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best integrated=skipped")
    if best_trace:
        b = best_trace
        print(
            f"  best trace={b['pairTransform']}/{b['spec']['kernel']}/"
            f"{b['spec']['center']}/{b['spec']['theta']}/lambda={b['spec']['lambda']} "
            f"{b['greenTraceGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best trace=skipped")
    print(f"  bridge constructed: {data['nonlocalDictionaryConstructedBridge']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
