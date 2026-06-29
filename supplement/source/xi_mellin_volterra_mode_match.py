#!/usr/bin/env python3
r"""Match finite-core Xi Mellin atoms to Volterra ratio atoms.

The generic coarea/resolvent probes try kernels that do not know the theta
mode structure.  This script tests the more rigid dictionary forced by the
finite-core expansion

    Phi(t) = sum_i a_i exp(lambda_i |t| - c_i exp(2|t|)).

For v=2t and Psi(v)=Phi(v/2), each positive-side atom is

    a_i exp(beta_i v - c_i exp(v)),       beta_i=lambda_i/2.

Its exact positive-side Xi contribution is

    (a_i/2) c_i^(-(beta_i+i z/2)) Gamma(beta_i+i z/2, c_i),

and the even Phi contribution is obtained by adding z -> -z.

The matching Volterra atom inside A_s(u)=Psi(s+u)/Psi(s) is

    rho_i(s) exp(beta_i u - c_i exp(s)(exp(u)-1)).

We use these paired atoms to build a spectral/mode-matched Hardy-to-Volterra
response and test whether it improves the joint Hardy branch Gram residual.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import default_z_nodes  # noqa: E402
from klm_debranges_continuous_normal_equation import (  # noqa: E402
    assemble_feature_forms,
    joint_from_forms,
    raw_lambda_traces,
    residual_report,
    trace_lift,
)
from klm_debranges_coupled_branch_normal_equation import PAIR_TRANSFORMS  # noqa: E402
from klm_debranges_trace_map_constructor import (  # noqa: E402
    ctranspose,
    finite_trace_model,
    frob_norm,
    hardy_targets,
    mapped_joint,
    to_float,
)
from klm_debranges_transmutation_kernel_probe import (  # noqa: E402
    lifted_feature_rows,
)
from mp_partial_k_corr import partial_phi, weights_for  # noqa: E402
from quotient_factorization_mp import endpoint_ratios  # noqa: E402
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402


def mode_weights(kind: str) -> dict[int, mp.mpf]:
    if kind == "raw1":
        return {1: mp.mpf("1")}
    return weights_for(kind)


def atom_records(kind: str) -> list[dict]:
    records = []
    for mode, weight in mode_weights(kind).items():
        c = mp.pi * mode * mode
        records.append(
            {
                "mode": mode,
                "term": "quartic",
                "coeff": weight * 4 * c * c,
                "beta": mp.mpf("2.25"),
                "c": c,
            }
        )
        records.append(
            {
                "mode": mode,
                "term": "quadratic",
                "coeff": weight * -6 * c,
                "beta": mp.mpf("1.25"),
                "c": c,
            }
        )
    return records


def positive_atom_xi(record: dict, z) -> mp.mpc:
    zz = mp.mpc(z)
    beta = mp.mpf(record["beta"])
    c = mp.mpf(record["c"])
    a = beta + mp.mpf("0.5") * mp.j * zz
    return mp.mpf("0.5") * record["coeff"] * c ** (-a) * mp.gammainc(a, c, mp.inf)


def atom_xi(record: dict, z, amp_kind: str) -> mp.mpc:
    if amp_kind == "positive":
        return positive_atom_xi(record, z)
    if amp_kind == "positive_reflected":
        return positive_atom_xi(record, -z)
    if amp_kind == "even":
        return positive_atom_xi(record, z) + positive_atom_xi(record, -z)
    raise ValueError(amp_kind)


class ExactFiniteCoreXi:
    def __init__(self, kind: str, amp_kind: str = "even"):
        self.kind = kind
        self.amp_kind = amp_kind
        self.records = atom_records(kind)

    def atom_values(self, z) -> list[mp.mpc]:
        return [atom_xi(record, z, self.amp_kind) for record in self.records]

    def __call__(self, z) -> mp.mpc:
        return mp.fsum(self.atom_values(z))


def simpson_grid_symmetric(tmax: mp.mpf, intervals: int):
    if intervals % 2:
        intervals += 1
    h = 2 * tmax / intervals
    pts = [-tmax + i * h for i in range(intervals + 1)]
    wts = []
    for i in range(intervals + 1):
        fac = 1 if i == 0 or i == intervals else 4 if i % 2 else 2
        wts.append(fac * h / 3)
    return pts, wts


def numeric_finite_core_xi(kind: str, z, tmax: mp.mpf, intervals: int) -> mp.mpc:
    weights = mode_weights(kind)
    pts, wts = simpson_grid_symmetric(tmax, intervals)
    zz = mp.mpc(z)
    return mp.fsum(w * partial_phi(t, weights) * mp.e ** (mp.j * zz * t) for t, w in zip(pts, wts))


def volterra_atom_values(s: mp.mpf, u: mp.mpf, pieces) -> list[mp.mpf]:
    es = mp.e**s
    eu = mp.e**u
    return [
        ratio * mp.e ** (beta * u - c * es * (eu - 1))
        for ratio, beta, c in endpoint_ratios(s, pieces)
    ]


def matched_response(args, xi_core: ExactFiniteCoreXi, z: complex, lifted: dict, branch: str, phase, target_kind: str, eps: int):
    omega = mp.mpf(args.omega)
    pieces = pieces_for(args.kind)
    z_shift = mp.mpc(z) + (mp.j * omega if branch == "plus" else -mp.j * omega)
    atom_amps = xi_core.atom_values(z_shift)
    y = mp.matrix(len(lifted["meta"]), 1)
    for i, (s, u, sigma, _A, root) in enumerate(lifted["meta"]):
        r = s + u
        branch_factor = mp.e ** (mp.mpf("0.5") * sigma * omega * r)
        atoms = volterra_atom_values(s, u, pieces)
        atom_sum = mp.fsum(atom_amps[j] * atoms[j] for j in range(len(atoms))) / mp.sqrt(2 * mp.pi)
        theta = mp.mpf("1")
        if target_kind == "atom":
            theta = mp.mpf("1")
        elif target_kind == "atom_branch":
            theta = branch_factor
        elif target_kind == "atom_feature":
            theta = branch_factor * (1 + eps * r)
        elif target_kind == "atom_coarea":
            theta = 1 / mp.sqrt(min(mp.mpf(args.L), max(r, mp.mpf("1e-30"))))
        elif target_kind == "atom_coarea_branch":
            theta = branch_factor / mp.sqrt(min(mp.mpf(args.L), max(r, mp.mpf("1e-30"))))
        elif target_kind == "atom_coarea_feature":
            theta = branch_factor * (1 + eps * r) / mp.sqrt(min(mp.mpf(args.L), max(r, mp.mpf("1e-30"))))
        else:
            raise ValueError(target_kind)
        y[i] = root * mp.mpc(phase) * atom_sum * theta
    return y


def solve_for_spec(args, xi_core: ExactFiniteCoreXi, lifted: dict, pair_transform: str, target_kind: str):
    nodes = default_z_nodes()
    A = ctranspose(lifted["plus"]) * lifted["plus"] + ctranspose(lifted["minus"]) * lifted["minus"]
    for i in range(A.rows):
        A[i, i] += mp.mpf(args.reg)
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    for j, z in enumerate(nodes):
        yp = matched_response(args, xi_core, z, lifted, plus_target[0], plus_target[1], target_kind, 1)
        ym = matched_response(args, xi_core, z, lifted, minus_target[0], minus_target[1], target_kind, -1)
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


def evaluate(args, xi_core, forms, lifted, trace_model, target, pair_transform: str, target_kind: str):
    X, residuals = solve_for_spec(args, xi_core, lifted, pair_transform, target_kind)
    lift = lifted_joint(lifted, X)
    row = {
        "pairTransform": pair_transform,
        "targetKind": target_kind,
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


def atom_xi_diagnostics(args, xi_core: ExactFiniteCoreXi):
    rows = []
    for z in default_z_nodes():
        exact = xi_core(z)
        numeric = numeric_finite_core_xi(args.kind, z, mp.mpf(args.xi_tmax), args.xi_intervals)
        rows.append(
            {
                "z": {"real": float(z.real), "imag": float(z.imag)},
                "exactReal": to_float(mp.re(exact)),
                "exactImag": to_float(mp.im(exact)),
                "numericReal": to_float(mp.re(numeric)),
                "numericImag": to_float(mp.im(numeric)),
                "relativeError": to_float(abs(exact - numeric) / max(abs(numeric), mp.mpf("1e-300"))),
            }
        )
    return rows


def parse_csv(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


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
    parser.add_argument("--xi-tmax", default="8")
    parser.add_argument("--xi-intervals", type=int, default=900)
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--skip-trace-lift", action="store_true")
    parser.add_argument("--skip-integrated-forms", action="store_true")
    parser.add_argument("--amp-kind", choices=["even", "positive", "positive_reflected"], default="even")
    parser.add_argument(
        "--target-kinds",
        default="atom,atom_branch,atom_feature,atom_coarea,atom_coarea_branch,atom_coarea_feature",
    )
    parser.add_argument("--pair-transforms", default="direct,swap")
    parser.add_argument("--json-out", default="xi_mellin_volterra_mode_match.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi_core = ExactFiniteCoreXi(args.kind, args.amp_kind)
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
    target = hardy_targets(xi_core, float(args.omega), None)
    pair_transforms = parse_csv(args.pair_transforms)
    target_kinds = parse_csv(args.target_kinds)
    results = []
    for pair in pair_transforms:
        if pair not in PAIR_TRANSFORMS:
            raise ValueError(pair)
        for target_kind in target_kinds:
            results.append(evaluate(args, xi_core, forms, lifted, trace_model, target, pair, target_kind))

    best_lifted = min(results, key=lambda row: row["liftedTriangleGram"]["scaledJointRelative"])
    integrated_results = [row for row in results if row["integratedVolterraGram"] is not None]
    best_integrated = (
        min(integrated_results, key=lambda row: row["integratedVolterraGram"]["scaledJointRelative"])
        if integrated_results
        else None
    )
    trace_results = [row for row in results if row["greenTraceGram"] is not None]
    best_trace = min(trace_results, key=lambda row: row["greenTraceGram"]["scaledJointRelative"]) if trace_results else None
    diagnostics = atom_xi_diagnostics(args, xi_core) if args.amp_kind == "even" else []
    max_xi_diag = max((row["relativeError"] for row in diagnostics), default=None)
    data = {
        "theoremName": "Xi Mellin atom to Volterra mode matching probe",
        "ansatz": "spectral atom split of Xi matched term-by-term to Volterra ratio atoms in A_s(u)",
        "kind": args.kind,
        "omega": float(args.omega),
        "basis": args.basis,
        "constraints": args.constraints,
        "ampKind": args.amp_kind,
        "atomCount": len(xi_core.records),
        "atomRecords": [
            {
                "mode": rec["mode"],
                "term": rec["term"],
                "coeff": to_float(rec["coeff"]),
                "beta": to_float(rec["beta"]),
                "c": to_float(rec["c"]),
            }
            for rec in xi_core.records
        ],
        "xiAtomDiagnostics": diagnostics,
        "maxXiAtomDiagnosticRelativeError": max_xi_diag,
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "integratedFormsSkipped": bool(args.skip_integrated_forms),
        "results": results,
        "bestLiftedTriangleGram": best_lifted,
        "bestIntegratedVolterraGram": best_integrated,
        "bestGreenTraceGram": best_trace,
        "modeMatchedBridgeConstructed": bool(
            best_trace and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6
        ),
        "diagnosis": (
            "The exact finite-core Xi atom expansion is matched to the exact "
            "Volterra ratio atoms.  If the Xi atom diagnostic is small but the "
            "Gram residual remains large, the missing transmutation is not a "
            "simple diagonal spectral mode split; it must also mix modes or "
            "carry another boundary/Hardy factor."
        ),
        "nextProofTarget": (
            "If mode matching does not close the joint Gram, derive the actual "
            "mode-mixing matrix between the Xi incomplete-gamma atoms and the "
            "Volterra ratio atoms, rather than a diagonal atom dictionary."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Xi Mellin atom to Volterra mode matching probe")
    print(
        f"  kind={args.kind} omega={args.omega} basis={args.basis} "
        f"constraints={args.constraints} atoms={len(xi_core.records)} amp={args.amp_kind}"
    )
    if max_xi_diag is not None:
        print(f"  max exact-atom Xi vs finite quadrature relative error: {max_xi_diag:.3e}")
    for row in results:
        integ = (
            f"{row['integratedVolterraGram']['scaledJointRelative']:.3e}"
            if row["integratedVolterraGram"] is not None
            else "skipped"
        )
        trace = (
            f"{row['greenTraceGram']['scaledJointRelative']:.3e}"
            if row["greenTraceGram"] is not None
            else "skipped"
        )
        print(
            f"  {row['pairTransform']:6s} {row['targetKind']:20s} "
            f"lifted={row['liftedTriangleGram']['scaledJointRelative']:.3e} "
            f"integrated={integ} trace={trace}"
        )
    b = best_lifted
    print(
        f"  best lifted={b['pairTransform']}/{b['targetKind']} "
        f"{b['liftedTriangleGram']['scaledJointRelative']:.6e}"
    )
    if best_integrated:
        b = best_integrated
        print(
            f"  best integrated={b['pairTransform']}/{b['targetKind']} "
            f"{b['integratedVolterraGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best integrated=skipped")
    if best_trace:
        b = best_trace
        print(
            f"  best trace={b['pairTransform']}/{b['targetKind']} "
            f"{b['greenTraceGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best trace=skipped")
    print(f"  bridge constructed: {data['modeMatchedBridgeConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
