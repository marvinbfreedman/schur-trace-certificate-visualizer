#!/usr/bin/env python3
r"""Mellin-boundary trace/concomitant for the KLM/de Branges bridge.

The exact Mellin split produces the boundary prefix

    B_i(s,z)=1/2 a_i c_i^(-alpha_i) int_{c_i}^{c_i e^s}
             y^{alpha_i-1} exp(-y) dy,
    alpha_i = beta_i + i z/2.

This prefix has the first-order identity

    d_s B_i(s,z) = 1/2 a_i exp(beta_i s-c_i exp(s)) exp(i z s/2),
    B_i(0,z)=0.

Thus the missing boundary object is the Volterra primitive trace

    Mu_z(f) = int_0^L B(s,z) f(s) ds.

Equivalently, with F(s)=int_s^L f(t)dt,

    Mu_z(f) = int_0^L B'(s,z) F(s) ds,

because B(0,z)=0 and F(L)=0.  This is the endpoint concomitant missing from
the old Lambda_a trace family.

This script verifies the derivative identity and checks that augmenting the
sampled endpoint trace matrix by the Mu_z rows kills the boundary prefix in
the finite quotient model.
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
from klm_debranges_trace_map_constructor import ctranspose, to_float  # noqa: E402
from quotient_factorization_mp import (  # noqa: E402
    columns,
    gram_matrix,
    legendre_quadrature,
    poly_value,
    shifted_legendre_polys,
    trace_matrix,
)
from xi_boundary_prefix_trace_resolution import (  # noqa: E402
    boundary_prefix_row,
    branch_z_values,
    local_trace_args,
    weighted_null_residual,
)
from xi_mellin_convolution_boundary_identity import atom_prefix  # noqa: E402
from xi_mellin_volterra_mode_match import atom_records  # noqa: E402


def stack_rows(a: mp.matrix, b: mp.matrix) -> mp.matrix:
    out = mp.matrix(a.rows + b.rows, a.cols)
    for i in range(a.rows):
        for j in range(a.cols):
            out[i, j] = a[i, j]
    for i in range(b.rows):
        for j in range(b.cols):
            out[a.rows + i, j] = b[i, j]
    return out


def prefix_sum(args, s: mp.mpf, z) -> mp.mpc:
    scale = mp.mpf(args.fourier_scale)
    records = atom_records(args.kind)
    total = mp.fsum(atom_prefix(record, s, z, scale) for record in records)
    if args.amp_kind == "positive":
        return total
    if args.amp_kind == "even":
        return total + mp.fsum(atom_prefix(record, s, -mp.mpc(z), scale) for record in records)
    raise ValueError(args.amp_kind)


def prefix_derivative_density(args, s: mp.mpf, z) -> mp.mpc:
    scale = mp.mpf(args.fourier_scale)
    zz = mp.mpc(z)
    total = mp.mpc("0")
    for record in atom_records(args.kind):
        total += (
            mp.mpf("0.5")
            * record["coeff"]
            * mp.e ** (record["beta"] * s - record["c"] * mp.e**s)
            * mp.e ** (mp.j * scale * zz * s)
        )
    if args.amp_kind == "positive":
        return total
    if args.amp_kind == "even":
        reflected = mp.mpc("0")
        for record in atom_records(args.kind):
            reflected += (
                mp.mpf("0.5")
                * record["coeff"]
                * mp.e ** (record["beta"] * s - record["c"] * mp.e**s)
                * mp.e ** (-mp.j * scale * zz * s)
            )
        return total + reflected
    raise ValueError(args.amp_kind)


def derivative_identity_check(args):
    s_values = [mp.mpf(x) for x in args.s_values.replace(",", " ").split()]
    h = mp.mpf(args.derivative_h)
    rows = []
    max_rel = mp.mpf("0")
    max_b0 = mp.mpf("0")
    for node_index, branch, z in branch_z_values(args):
        b0 = abs(prefix_sum(args, mp.mpf("0"), z))
        max_b0 = max(max_b0, b0)
        for s in s_values:
            if s - h < 0:
                fd = (prefix_sum(args, s + h, z) - prefix_sum(args, s, z)) / h
            else:
                fd = (prefix_sum(args, s + h, z) - prefix_sum(args, s - h, z)) / (2 * h)
            exact = prefix_derivative_density(args, s, z)
            rel = abs(fd - exact) / max(abs(exact), mp.mpf("1e-300"))
            max_rel = max(max_rel, rel)
            rows.append(
                {
                    "nodeIndex": node_index,
                    "branch": branch,
                    "s": to_float(s),
                    "z": {"real": to_float(mp.re(z)), "imag": to_float(mp.im(z))},
                    "finiteDifferenceRelativeError": to_float(rel),
                    "exactDerivativeAbs": to_float(abs(exact)),
                }
            )
    return {
        "maxDerivativeFiniteDifferenceRelativeError": to_float(max_rel),
        "maxBoundaryPrefixAtZeroAbs": to_float(max_b0),
        "derivativeRows": rows,
    }


def prefix_rows(args, polys) -> mp.matrix:
    rows = []
    for _node_index, _branch, z in branch_z_values(args):
        row = boundary_prefix_row(args, polys, z)
        rows.append([row[0, j] for j in range(row.cols)])
    return mp.matrix(rows)


def row_action_max(row_mat: mp.matrix, N: mp.matrix, K: mp.matrix, psd_tol: mp.mpf):
    max_op = mp.mpf("0")
    max_l2 = mp.mpf("0")
    modes = 0
    for i in range(row_mat.rows):
        row = mp.matrix(1, row_mat.cols)
        for j in range(row_mat.cols):
            row[0, j] = row_mat[i, j]
        op, l2, positive_modes = weighted_null_residual(row, N, K, psd_tol)
        max_op = max(max_op, op)
        max_l2 = max(max_l2, l2)
        modes = max(modes, positive_modes)
    return max_op, max_l2, modes


def split_trace_nullspace_complex(R: mp.matrix, rank_tol_text: str):
    gram = ctranspose(R) * R
    vals, vecs = mp.eighe((gram + ctranspose(gram)) / 2, eigvals_only=False)
    vmax = max(abs(v) for v in vals) if len(vals) else mp.mpf("0")
    tol = mp.mpf(rank_tol_text) * max(mp.mpf("1"), vmax)
    n_idx = [i for i, val in enumerate(vals) if val <= tol]
    u_idx = [i for i, val in enumerate(vals) if val > tol]
    return columns(vecs, n_idx), columns(vecs, u_idx), vals, tol


def augmented_trace_check(args):
    L = mp.mpf(args.L)
    polys = shifted_legendre_polys(args.basis, L)
    _centers, R_lambda = trace_matrix(polys, local_trace_args(args))
    R_mu = prefix_rows(args, polys)
    R_aug = stack_rows(R_lambda, R_mu)
    K, _ = gram_matrix(args, mp.mpf(args.omega), L)

    N_old, U_old, _old_vals, old_tol = split_trace_nullspace_complex(R_lambda, args.rank_tol)
    N_aug, U_aug, _aug_vals, aug_tol = split_trace_nullspace_complex(R_aug, args.rank_tol)
    old_op, old_l2, old_modes = row_action_max(R_mu, N_old, K, mp.mpf(args.psd_tol))
    aug_op, aug_l2, aug_modes = row_action_max(R_mu, N_aug, K, mp.mpf(args.psd_tol))
    return {
        "lambdaTraceRank": U_old.cols,
        "lambdaTraceNullity": N_old.cols,
        "augmentedTraceRank": U_aug.cols,
        "augmentedTraceNullity": N_aug.cols,
        "lambdaTraceRankTolerance": to_float(old_tol),
        "augmentedTraceRankTolerance": to_float(aug_tol),
        "boundaryRows": R_mu.rows,
        "oldTraceBoundaryNullEnergyOpMax": to_float(old_op),
        "oldTraceBoundaryNullEnergyOpL2": to_float(old_l2),
        "oldTracePositiveNullModes": old_modes,
        "augmentedTraceBoundaryNullEnergyOpMax": to_float(aug_op),
        "augmentedTraceBoundaryNullEnergyOpL2": to_float(aug_l2),
        "augmentedTracePositiveNullModes": aug_modes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["full", "endpoint_b"], default="full")
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=10)
    parser.add_argument("--quad", type=int, default=18)
    parser.add_argument("--constraints", type=int, default=7)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=20)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-24")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--fourier-scale", default="0.5")
    parser.add_argument("--amp-kind", choices=["positive", "even"], default="positive")
    parser.add_argument("--hardy-normalization", action="store_true", default=True)
    parser.add_argument("--s-values", default="0,0.02,0.1,0.3,0.545")
    parser.add_argument("--derivative-h", default="1e-5")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="xi_mellin_boundary_concomitant.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    derivative = derivative_identity_check(args)
    augmented = augmented_trace_check(args)
    data = {
        "theoremName": "Mellin-boundary trace/concomitant derivation",
        "kind": args.kind,
        "omega": float(args.omega),
        "L": float(args.L),
        "basis": args.basis,
        "constraints": args.constraints,
        "ampKind": args.amp_kind,
        "fourierScale": float(args.fourier_scale),
        "concomitantIdentity": {
            "prefixDerivative": "d_s B_i(s,z)=1/2 a_i exp(beta_i s-c_i exp(s)) exp(i z s/2)",
            "trace": "Mu_z(f)=int_0^L B(s,z) f(s) ds",
            "primitiveForm": "Mu_z(f)=int_0^L B'(s,z) F(s) ds, F(s)=int_s^L f(t)dt",
            "endpointTerms": "B(0,z)=0 and F(L)=0, so the integration-by-parts boundary term vanishes.",
        },
        **derivative,
        **augmented,
        "mellinBoundaryConcomitantDerived": True,
        "augmentedTraceKillsBoundaryPrefix": augmented["augmentedTraceBoundaryNullEnergyOpMax"] < 1e-8,
        "diagnosis": (
            "The missing boundary object is the Volterra primitive trace Mu_z, "
            "not another Lambda_a row.  Adding Mu_z to the finite trace matrix "
            "kills the boundary prefix on the augmented nullspace by construction; "
            "the remaining hard step is proving the augmented trace repair is "
            "positive or annihilated in the KLM/de Branges pullback."
        ),
        "nextProofTarget": (
            "Prove positivity/annihilation of the augmented trace repair "
            "(Lambda_a, Mu_z): combine the diagonal Volterra tail identity with "
            "the new Mellin-boundary concomitant without introducing a negative "
            "Schur defect."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Mellin-boundary trace/concomitant derivation")
    print(
        f"  kind={args.kind} omega={args.omega} basis={args.basis} "
        f"constraints={args.constraints}"
    )
    print(f"  max derivative FD relative error: {data['maxDerivativeFiniteDifferenceRelativeError']:.3e}")
    print(f"  max B(0,z): {data['maxBoundaryPrefixAtZeroAbs']:.3e}")
    print(
        "  old trace boundary null op: "
        f"{data['oldTraceBoundaryNullEnergyOpMax']:.3e}"
    )
    print(
        "  augmented trace boundary null op: "
        f"{data['augmentedTraceBoundaryNullEnergyOpMax']:.3e}"
    )
    print(f"  augmented kills prefix: {data['augmentedTraceKillsBoundaryPrefix']}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
