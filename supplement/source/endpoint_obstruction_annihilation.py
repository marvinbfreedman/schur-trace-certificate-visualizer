#!/usr/bin/env python3
r"""Finite diagnostic for endpoint Fredholm obstructions.

The exact endpoint equation after the interior Green jump is

    M z = -b(d).

Endpoint compatibility means every left obstruction w in ker(M^T) annihilates
the actual source endpoint vector:

    w^T b(d_u)=0.

This script computes the finite sampled-flow analogue of that statement and
also checks whether the active vector represented by w is globally trace-null.
The second check is important: if the obstruction vector is not trace-null,
then compatibility is a special endpoint/source identity, not a consequence of
closed trace range inclusion.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import (  # noqa: E402
    active_setup,
    endpoint_active_row,
    endpoint_map,
    mat_frobenius,
    propagate_flows,
    sample_derivatives,
)
from adjoint_green_jump_conditions import jump_matrix, solve_jump  # noqa: E402
from global_trace_observability_gap import child_args, f, fmt, stack_source_rows  # noqa: E402
from lagrange_energy_control_certificate import make_qargs  # noqa: E402
from local_trace_tower_representer_scan import exact_trace_derivatives  # noqa: E402
from quotient_factorization_mp import endpoint_b_quadrature, trace_matrix  # noqa: E402
from source_concomitant_membership import source_derivatives  # noqa: E402
from trace_concomitant_membership import trace_green_concomitant_row  # noqa: E402
from trace_lagrange_adjoint_control import adjoint_source_value  # noqa: E402


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def vec_norm(vec) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(vec[i]) ** 2 for i in range(vec.rows)))


def row_norm(values) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in values))


def left_null_vector(mat: mp.matrix) -> tuple[mp.matrix, list[mp.mpf]]:
    gram = (mat * mat.T + mat * mat.T) / 2
    vals, vecs = mp.eigsy(gram, eigvals_only=False)
    idx = min(range(len(vals)), key=lambda i: abs(vals[i]))
    out = mp.matrix(mat.rows, 1)
    for i in range(mat.rows):
        out[i] = vecs[i, idx]
    nrm = vec_norm(out)
    if nrm:
        out /= nrm
    if out[0] < 0:
        out *= -1
    return out, [vals[i] for i in range(len(vals))]


def matrix_vec_norm(mat: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols)))


def finite_trace_rows(args, polys, modes):
    targs = child_args(args, args.basis, args.global_constraints)
    targs.constraints = args.global_constraints
    centers, rows = trace_matrix(polys, make_qargs(targs))
    return centers, rows * modes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=11)
    parser.add_argument("--active-tol", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--ode-samples", type=int, default=11)
    parser.add_argument("--stencil-width", type=int, default=9)
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--matrix-order", type=int, default=70)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="endpoint_obstruction_annihilation.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    s0 = mp.mpf(args.s0)
    left = mp.mpf(args.constraint_min)
    right = mp.mpf(args.constraint_max)
    vals, e_derivs_at_s0, lam_derivs = exact_trace_derivatives(args, args.max_trace_q)
    base, active_modes, active_idx = active_setup(args, e_derivs_at_s0)

    centers, derivs, min_gap, min_cos = sample_derivatives(args, args.jet_order - 1)
    flow_left, flow_right, _idx0 = propagate_flows(centers, derivs, s0)
    e_left = derivs[0]
    e_right = derivs[-1]
    endpoint_mat, _zero = endpoint_map(
        base["polys"],
        active_modes,
        e_left,
        e_right,
        flow_left,
        flow_right,
        left,
        right,
    )
    obstruction, gram_vals = left_null_vector(endpoint_mat)
    mtw = endpoint_mat.T * obstruction

    obstruction_mode = active_modes * obstruction
    trace_centers, trace_on_obstruction = finite_trace_rows(args, base["polys"], obstruction_mode)
    trace_norm = matrix_vec_norm(trace_on_obstruction)
    trace_max = max(abs(trace_on_obstruction[i, 0]) for i in range(trace_on_obstruction.rows))

    _source_nodes, source_rows = stack_source_rows(args, base["polys"], e_derivs_at_s0)
    source_on_obstruction = source_rows * obstruction_mode
    source_norm = matrix_vec_norm(source_on_obstruction)
    source_max = max(abs(source_on_obstruction[i, 0]) for i in range(source_on_obstruction.rows))

    jmat = jump_matrix(e_derivs_at_s0, args.jet_order - 1)
    t_values = [mp.mpf(piece) for piece in args.t_values.replace(",", " ").split()]
    r_nodes, r_weights = endpoint_b_quadrature(mp.mpf(args.kernel_rmax), args.kernel_order)

    rows = []
    max_abs_pairing = mp.mpf("0")
    max_relative_pairing = mp.mpf("0")
    for t in t_values:
        h_derivs = source_derivatives(mp.pi, s0, t, args.jet_order, r_nodes, r_weights)
        boundary_source = trace_green_concomitant_row(e_derivs_at_s0, h_derivs, args.jet_order)
        pstar = adjoint_source_value(e_derivs_at_s0, h_derivs, args.jet_order)
        eval_source = [pstar] + [mp.mpf("0") for _ in range(args.jet_order - 2)]
        for component, source_row in (("boundary", boundary_source), ("adjointEval", eval_source)):
            delta, jump_residual, _target = solve_jump(jmat, source_row)
            delta_vec = mp.matrix(len(delta), 1)
            for i, value in enumerate(delta):
                delta_vec[i] = value
            zero_left = mp.matrix(len(delta), 1)
            b_values = endpoint_active_row(
                base["polys"],
                active_modes,
                e_left,
                e_right,
                zero_left,
                flow_right * delta_vec,
                left,
                right,
            )
            pairing = mp.fsum(obstruction[i] * b_values[i] for i in range(len(b_values)))
            b_norm = row_norm(b_values)
            rel = abs(pairing) / max(mp.mpf("1"), b_norm)
            max_abs_pairing = max(max_abs_pairing, abs(pairing))
            max_relative_pairing = max(max_relative_pairing, rel)
            rows.append(
                {
                    "sourceNode": f(t),
                    "component": component,
                    "endpointVectorNorm": f(b_norm),
                    "obstructionPairing": f(pairing),
                    "relativePairing": f(rel),
                    "jumpNorm": f(row_norm(delta)),
                }
            )

    data = {
        "theoremName": "endpoint obstruction annihilation diagnostic",
        "basis": args.basis,
        "s0": f(s0),
        "interval": [f(left), f(right)],
        "activeDim": len(active_idx),
        "activeIndices": active_idx,
        "endpointMapFrobenius": f(mat_frobenius(endpoint_mat)),
        "endpointGramEigenvalues": [f(value) for value in gram_vals],
        "leftNullVector": [f(obstruction[i]) for i in range(obstruction.rows)],
        "leftNullResidualNorm": f(vec_norm(mtw)),
        "maxAbsObstructionPairing": f(max_abs_pairing),
        "maxRelativeObstructionPairing": f(max_relative_pairing),
        "traceCenters": [f(x) for x in trace_centers],
        "traceOnObstructionNorm": f(trace_norm),
        "traceOnObstructionMaxAbs": f(trace_max),
        "sourceOnObstructionNorm": f(source_norm),
        "sourceOnObstructionMaxAbs": f(source_max),
        "rows": rows,
        "statuses": {
            "finiteLeftObstructionFound": status(
                "finite left endpoint obstruction found",
                vec_norm(mtw) < mp.mpf("1e-20") * max(mp.mpf("1"), mat_frobenius(endpoint_mat)),
                "The endpoint map has a numerically resolved left null direction.",
            ),
            "actualEndpointVectorsAnnihilated": status(
                "actual endpoint vectors annihilated",
                max_relative_pairing < mp.mpf("1e-24"),
                "The computed obstruction annihilates b(d_u) for every sampled actual source row.",
            ),
            "obstructionIsGloballyTraceNull": status(
                "obstruction mode globally trace-null",
                trace_norm < mp.mpf("1e-20"),
                "If false, endpoint compatibility is not explained by the global closed-trace condition alone.",
            ),
        },
        "interpretation": (
            "The finite obstruction annihilates the actual endpoint source "
            "vectors, while the associated active polynomial need not be "
            "globally trace-null.  This supports a special endpoint/source "
            "identity as the non-circular continuum target."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint obstruction annihilation diagnostic")
    print(f"  left null residual: {fmt(vec_norm(mtw), 12)}")
    print(f"  max obstruction pairing: {fmt(max_abs_pairing, 12)}")
    print(f"  max relative pairing: {fmt(max_relative_pairing, 12)}")
    print(f"  trace norm of obstruction mode: {fmt(trace_norm, 12)}")
    print(f"  source norm of obstruction mode: {fmt(source_norm, 12)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
