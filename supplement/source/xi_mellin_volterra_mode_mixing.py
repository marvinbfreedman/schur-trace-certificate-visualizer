#!/usr/bin/env python3
r"""Finite mode-mixing matrix between Xi atoms and Volterra atoms.

The diagonal atom dictionary

    Xi atom i -> Volterra ratio atom i

is normalized correctly but fails the joint Hardy/Volterra Gram test.  This
script derives the next finite object: a 6x6 matrix M such that, on the lifted
Volterra triangle,

    B_eps(s,u,sigma) M a_branch(z) ~= h_branch,z(s,u,sigma).

Here a_branch(z) is the vector of incomplete-gamma Xi atom amplitudes and
B_eps is the row matrix of the six Volterra ratio atoms.  The matrix M is
learned from the row-level least-squares identity, then tested against the
joint Hardy branch Gram without refitting the Gram itself.
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
from klm_debranges_transmutation_kernel_probe import lifted_feature_rows  # noqa: E402
from quotient_factorization_mp import endpoint_ratios  # noqa: E402
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402
from xi_mellin_volterra_mode_match import ExactFiniteCoreXi, atom_xi_diagnostics  # noqa: E402


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


def theta_value(kind: str, args, s: mp.mpf, u: mp.mpf, sigma: int, eps: int) -> mp.mpf:
    r = s + u
    omega = mp.mpf(args.omega)
    branch = mp.e ** (mp.mpf("0.5") * sigma * omega * r)
    coarea = 1 / mp.sqrt(min(mp.mpf(args.L), max(r, mp.mpf("1e-30"))))
    if kind == "plain":
        return mp.mpf("1")
    if kind == "branch":
        return branch
    if kind == "feature":
        return branch * (1 + eps * r)
    if kind == "coarea":
        return coarea
    if kind == "coarea_branch":
        return coarea * branch
    if kind == "coarea_feature":
        return coarea * branch * (1 + eps * r)
    raise ValueError(kind)


def volterra_atom_values(s: mp.mpf, u: mp.mpf, pieces) -> list[mp.mpf]:
    es = mp.e**s
    eu = mp.e**u
    return [
        ratio * mp.e ** (beta * u - c * es * (eu - 1))
        for ratio, beta, c in endpoint_ratios(s, pieces)
    ]


def atom_row_matrix(args, lifted: dict, theta_kind: str, eps: int) -> mp.matrix:
    pieces = pieces_for(args.kind)
    rows = mp.matrix(len(lifted["meta"]), len(pieces))
    for i, (s, u, sigma, _A, root) in enumerate(lifted["meta"]):
        theta = theta_value(theta_kind, args, s, u, sigma, eps)
        atoms = volterra_atom_values(s, u, pieces)
        for j, atom in enumerate(atoms):
            rows[i, j] = root * theta * atom
    return rows


def atom_amplitudes(args, xi_core: ExactFiniteCoreXi, z: complex, branch: str, phase) -> mp.matrix:
    omega = mp.mpf(args.omega)
    z_shift = mp.mpc(z) + (mp.j * omega if branch == "plus" else -mp.j * omega)
    vals = xi_core.atom_values(z_shift)
    out = mp.matrix(len(vals), 1)
    for i, val in enumerate(vals):
        out[i] = mp.mpc(phase) * val / mp.sqrt(2 * mp.pi)
    return out


def hardy_target_vector(args, xi_core: ExactFiniteCoreXi, z: complex, lifted: dict, branch: str, phase, theta_kind: str, center: str, eps: int):
    omega = mp.mpf(args.omega)
    z_shift = mp.mpc(z) + (mp.j * omega if branch == "plus" else -mp.j * omega)
    amp = mp.mpc(phase) * xi_core(z_shift) / mp.sqrt(2 * mp.pi)
    out = mp.matrix(len(lifted["meta"]), 1)
    for i, (s, u, sigma, _A, root) in enumerate(lifted["meta"]):
        c = center_value(center, s, u)
        theta = theta_value(theta_kind, args, s, u, sigma, eps)
        out[i] = root * theta * amp * mp.e ** (mp.j * mp.mpc(z) * c)
    return out


def vec_index(out_atom: int, in_atom: int, atom_count: int) -> int:
    return out_atom * atom_count + in_atom


def build_mixing_system(args, xi_core: ExactFiniteCoreXi, lifted: dict, pair_transform: str, theta_kind: str, center: str):
    nodes = default_z_nodes()
    atom_count = len(xi_core.records)
    rows_plus = atom_row_matrix(args, lifted, theta_kind, 1)
    rows_minus = atom_row_matrix(args, lifted, theta_kind, -1)
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    total_rows = 2 * len(nodes) * len(lifted["meta"])
    K = mp.matrix(total_rows, atom_count * atom_count)
    y = mp.matrix(total_rows, 1)
    row0 = 0
    samples = [
        (rows_plus, 1, plus_target[0], plus_target[1]),
        (rows_minus, -1, minus_target[0], minus_target[1]),
    ]
    for z in nodes:
        for B, eps, branch, phase in samples:
            a = atom_amplitudes(args, xi_core, z, branch, phase)
            target = hardy_target_vector(args, xi_core, z, lifted, branch, phase, theta_kind, center, eps)
            for r in range(B.rows):
                y[row0 + r] = target[r]
                for out_atom in range(atom_count):
                    brow = B[r, out_atom]
                    if brow == 0:
                        continue
                    for in_atom in range(atom_count):
                        K[row0 + r, vec_index(out_atom, in_atom, atom_count)] = brow * a[in_atom]
            row0 += B.rows
    return K, y, rows_plus, rows_minus


def solve_mixing_matrix(args, xi_core: ExactFiniteCoreXi, lifted: dict, pair_transform: str, theta_kind: str, center: str):
    K, y, rows_plus, rows_minus = build_mixing_system(args, xi_core, lifted, pair_transform, theta_kind, center)
    H = ctranspose(K) * K
    g = ctranspose(K) * y
    diag_max = max([abs(H[i, i]) for i in range(H.rows)] + [mp.mpf("1")])
    reg = mp.mpf(args.mix_reg) * diag_max
    for i in range(H.rows):
        H[i, i] += reg
    x = mp.lu_solve(H, g)
    pred = K * x
    residual = frob_norm(pred - y) / max(frob_norm(y), mp.mpf("1e-300"))
    atom_count = len(xi_core.records)
    M = mp.matrix(atom_count, atom_count)
    for out_atom in range(atom_count):
        for in_atom in range(atom_count):
            M[out_atom, in_atom] = x[vec_index(out_atom, in_atom, atom_count)]
    herm_H = (H + ctranspose(H)) / 2
    try:
        gram_vals = mp.eighe(herm_H, eigvals_only=True)
    except Exception:
        gram_vals = []
    cond = abs(gram_vals[-1]) / max(abs(gram_vals[0]), mp.mpf("1e-300")) if len(gram_vals) else mp.inf
    return {
        "M": M,
        "rowsPlus": rows_plus,
        "rowsMinus": rows_minus,
        "fitResidual": to_float(residual),
        "normalConditionEstimate": to_float(cond),
        "normalEigenMin": to_float(gram_vals[0]) if len(gram_vals) else None,
        "normalEigenMax": to_float(gram_vals[-1]) if len(gram_vals) else None,
    }


def mixed_response(args, xi_core: ExactFiniteCoreXi, z: complex, rows: mp.matrix, M: mp.matrix, branch: str, phase):
    a = atom_amplitudes(args, xi_core, z, branch, phase)
    return rows * M * a


def solve_coefficients_for_mixing(args, xi_core, lifted, mix: dict, pair_transform: str):
    nodes = default_z_nodes()
    A = ctranspose(lifted["plus"]) * lifted["plus"] + ctranspose(lifted["minus"]) * lifted["minus"]
    for i in range(A.rows):
        A[i, i] += mp.mpf(args.reg)
    X = mp.matrix(args.basis, len(nodes))
    residuals = []
    plus_target, minus_target = PAIR_TRANSFORMS[pair_transform]
    for j, z in enumerate(nodes):
        yp = mixed_response(args, xi_core, z, mix["rowsPlus"], mix["M"], plus_target[0], plus_target[1])
        ym = mixed_response(args, xi_core, z, mix["rowsMinus"], mix["M"], minus_target[0], minus_target[1])
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


def matrix_to_json(M: mp.matrix):
    rows = []
    for i in range(M.rows):
        row = []
        for j in range(M.cols):
            row.append({"real": to_float(mp.re(M[i, j])), "imag": to_float(mp.im(M[i, j]))})
        rows.append(row)
    return rows


def evaluate_spec(args, xi_core, forms, lifted, trace_model, target, pair_transform: str, theta_kind: str, center: str):
    mix = solve_mixing_matrix(args, xi_core, lifted, pair_transform, theta_kind, center)
    X, normal_residuals = solve_coefficients_for_mixing(args, xi_core, lifted, mix, pair_transform)
    lift = lifted_joint(lifted, X)
    row = {
        "pairTransform": pair_transform,
        "theta": theta_kind,
        "center": center,
        "mixingFitResidual": mix["fitResidual"],
        "mixingNormalConditionEstimate": mix["normalConditionEstimate"],
        "mixingNormalEigenMin": mix["normalEigenMin"],
        "mixingNormalEigenMax": mix["normalEigenMax"],
        "normalResidualMax": max(normal_residuals),
        "normalResidualMean": sum(normal_residuals) / len(normal_residuals),
        "coefficientNormFrobenius": to_float(frob_norm(X)),
        "mixingMatrixFrobenius": to_float(frob_norm(mix["M"])),
        "mixingMatrix": matrix_to_json(mix["M"]),
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
    parser.add_argument("--s-order", type=int, default=8)
    parser.add_argument("--u-order", type=int, default=16)
    parser.add_argument("--u-max", default="10")
    parser.add_argument("--laguerre-order", type=int, default=24)
    parser.add_argument("--rhs-laguerre-order", type=int, default=32)
    parser.add_argument("--xi-tmax", default="8")
    parser.add_argument("--xi-intervals", type=int, default=500)
    parser.add_argument("--reg", default="1e-30")
    parser.add_argument("--mix-reg", default="1e-28")
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--skip-trace-lift", action="store_true")
    parser.add_argument("--skip-integrated-forms", action="store_true")
    parser.add_argument("--amp-kind", choices=["even", "positive", "positive_reflected"], default="even")
    parser.add_argument("--thetas", default="plain,coarea,branch,coarea_branch")
    parser.add_argument("--centers", default="su,mellin_eta,u")
    parser.add_argument("--pair-transforms", default="direct,swap")
    parser.add_argument("--json-out", default="xi_mellin_volterra_mode_mixing.json")
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
    results = []
    for pair in parse_csv(args.pair_transforms):
        if pair not in PAIR_TRANSFORMS:
            raise ValueError(pair)
        for theta in parse_csv(args.thetas):
            for center in parse_csv(args.centers):
                results.append(evaluate_spec(args, xi_core, forms, lifted, trace_model, target, pair, theta, center))

    best_fit = min(results, key=lambda row: row["mixingFitResidual"])
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
        "theoremName": "finite Xi-to-Volterra atom mode-mixing matrix probe",
        "ansatz": "least-squares 6x6 matrix M from incomplete-gamma Xi atom amplitudes to Volterra ratio atoms",
        "kind": args.kind,
        "omega": float(args.omega),
        "basis": args.basis,
        "constraints": args.constraints,
        "ampKind": args.amp_kind,
        "atomCount": len(xi_core.records),
        "maxXiAtomDiagnosticRelativeError": max_xi_diag,
        "traceRank": trace_model["traceRank"] if trace_model is not None else None,
        "traceLiftSkipped": bool(args.skip_trace_lift),
        "integratedFormsSkipped": bool(args.skip_integrated_forms),
        "results": results,
        "bestMixingFit": best_fit,
        "bestLiftedTriangleGram": best_lifted,
        "bestIntegratedVolterraGram": best_integrated,
        "bestGreenTraceGram": best_trace,
        "modeMixingBridgeConstructed": bool(
            best_trace and best_trace["greenTraceGram"]["scaledJointRelative"] < 1e-6
        ),
        "diagnosis": (
            "This derives the finite 6x6 atom mixing matrix by row-level "
            "least squares, then tests the induced branch Gram.  A small row "
            "fit but large Gram residual means the chosen row target is not "
            "the canonical Hardy-to-Volterra transmutation; a large row fit "
            "means the six Volterra atoms cannot represent the chosen Hardy "
            "row target uniformly."
        ),
        "nextProofTarget": (
            "If this finite mode-mixing matrix does not close the bridge, derive "
            "the mixing from the exact Mellin convolution identity itself, not "
            "from row-level least squares on a coarea target."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Finite Xi-to-Volterra atom mode-mixing matrix probe")
    print(
        f"  kind={args.kind} omega={args.omega} basis={args.basis} "
        f"constraints={args.constraints} atoms={len(xi_core.records)}"
    )
    if max_xi_diag is not None:
        print(f"  max exact-atom Xi diagnostic relative error: {max_xi_diag:.3e}")
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
            f"  {row['pairTransform']:6s} {row['theta']:14s} {row['center']:10s} "
            f"fit={row['mixingFitResidual']:.3e} "
            f"lifted={row['liftedTriangleGram']['scaledJointRelative']:.3e} "
            f"integrated={integ} trace={trace}"
        )
    b = best_fit
    print(
        f"  best fit={b['pairTransform']}/{b['theta']}/{b['center']} "
        f"{b['mixingFitResidual']:.6e}"
    )
    b = best_lifted
    print(
        f"  best lifted={b['pairTransform']}/{b['theta']}/{b['center']} "
        f"{b['liftedTriangleGram']['scaledJointRelative']:.6e}"
    )
    if best_integrated:
        b = best_integrated
        print(
            f"  best integrated={b['pairTransform']}/{b['theta']}/{b['center']} "
            f"{b['integratedVolterraGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best integrated=skipped")
    if best_trace:
        b = best_trace
        print(
            f"  best trace={b['pairTransform']}/{b['theta']}/{b['center']} "
            f"{b['greenTraceGram']['scaledJointRelative']:.6e}"
        )
    else:
        print("  best trace=skipped")
    print(f"  bridge constructed: {data['modeMixingBridgeConstructed']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
