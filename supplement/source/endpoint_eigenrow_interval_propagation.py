#!/usr/bin/env python3
r"""Interval propagation through the confluent eigenrow Taylor recurrence.

The endpoint Krawczyk proof uses derivatives of the moving negative eigenrow
``e(s)`` of the endpoint confluent matrix ``A(s)``.  The previous hardening
step bounded the quadrature error for the Taylor matrices ``A_p``.  This script
propagates that matrix interval through the finite eigenvalue/eigenvector
Taylor system

    sum_{p=0}^n A_p v_{n-p} = sum_{p=0}^n lambda_p v_{n-p},
    sum_{p=0}^n <v_p,v_{n-p}> = delta_{n0},        0 <= n <= m.

It builds the exact Jacobian of this polynomial system at the computed center
and applies a small Krawczyk/Newton ball test.  This replaces the deliberately
huge placeholder amplification used by
``endpoint_confluent_segment_bernstein_certificate.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import status  # noqa: E402
from endpoint_confluent_trace_tail_certificate import collocation_nodes  # noqa: E402
from endpoint_kb_confluent_mp import integrate  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402
from trace_concomitant_exact_derivatives import center_taylor_mats, eigen_taylor  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def node_key(value: mp.mpf) -> str:
    return mp.nstr(value, 80)


def mat_to_text(mat: mp.matrix) -> list[list[str]]:
    return [[mp.nstr(mat[i, j], 90) for j in range(mat.cols)] for i in range(mat.rows)]


def text_to_mat(rows: list[list[str]]) -> mp.matrix:
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def vec_inf(values) -> mp.mpf:
    return max(abs(value) for value in values) if values else mp.mpf("0")


def mat_inf_norm(mat: mp.matrix) -> mp.mpf:
    return max(
        mp.fsum(abs(mat[i, j]) for j in range(mat.cols))
        for i in range(mat.rows)
    )


def flatten_index(n: int, kind: str, i: int, dim: int) -> int:
    base = n * (dim + 1)
    return base if kind == "lambda" else base + 1 + i


def build_jacobian(a_coeffs: list[mp.matrix], vectors, lambdas) -> mp.matrix:
    max_q = len(a_coeffs) - 1
    dim = a_coeffs[0].rows
    size = (max_q + 1) * (dim + 1)
    jac = mp.matrix(size, size)
    row = 0
    for n in range(max_q + 1):
        for i in range(dim):
            for k in range(n + 1):
                block = a_coeffs[n - k].copy()
                for ii in range(dim):
                    block[ii, ii] -= lambdas[n - k]
                for j in range(dim):
                    jac[row, flatten_index(k, "vector", j, dim)] += block[i, j]
            for k in range(n + 1):
                for_lambda = -vectors[n - k][i]
                jac[row, flatten_index(k, "lambda", 0, dim)] += for_lambda
            row += 1
        for k in range(n + 1):
            other = n - k
            for j in range(dim):
                jac[row, flatten_index(k, "vector", j, dim)] += vectors[other][j]
        row += 1
    return jac


def build_center(args, s: mp.mpf, max_q: int):
    cache_path = None
    if args.cache_dir:
        cache_dir = Path(args.cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = {
            "s": node_key(s),
            "dps": args.dps,
            "matrixOrder": args.condition_matrix_order,
            "matrixRmax": args.matrix_rmax,
            "jetOrder": args.jet_order,
            "maxQ": max_q,
            "kind": args.kind,
        }
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:24]
        cache_path = cache_dir / f"eigenrow_center_{digest}.json"
        if cache_path.exists():
            cached = load_json(str(cache_path))
            mats = [text_to_mat(rows) for rows in cached["aCoefficients"]]
            vals = [mp.mpf(value) for value in cached["eigenvalues"]]
            e_derivs = [
                [mp.mpf(value) for value in row]
                for row in cached["eDerivs"]
            ]
            lam_derivs = [mp.mpf(value) for value in cached["lambdaDerivs"]]
            return mats, vals, e_derivs, lam_derivs, True

    big_order = args.jet_order + max_q
    big, _segments = integrate(
        args.kind,
        mp.pi,
        s,
        big_order,
        mp.mpf(args.matrix_rmax),
        args.condition_matrix_order,
    )
    mats = center_taylor_mats(big, args.jet_order, max_q)
    vals, e_derivs, lam_derivs = eigen_taylor(mats, max_q)
    if cache_path is not None:
        cache_path.write_text(
            json.dumps(
                {
                    "s": node_key(s),
                    "aCoefficients": [mat_to_text(mat) for mat in mats],
                    "eigenvalues": [mp.nstr(value, 90) for value in vals],
                    "eDerivs": [
                        [mp.nstr(value, 90) for value in row]
                        for row in e_derivs
                    ],
                    "lambdaDerivs": [mp.nstr(value, 90) for value in lam_derivs],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return mats, vals, e_derivs, lam_derivs, False


def residual_inf_bound(vectors, entry_epsilon: mp.mpf) -> mp.mpf:
    max_q = len(vectors) - 1
    out = mp.mpf("0")
    for n in range(max_q + 1):
        row_residual = entry_epsilon * mp.fsum(
            mp.fsum(abs(value) for value in vectors[n - p])
            for p in range(n + 1)
        )
        out = max(out, row_residual)
    return out


def derivative_variation_bound(
    solution_radius: mp.mpf,
    entry_epsilon: mp.mpf,
    dim: int,
    max_q: int,
) -> mp.mpf:
    eig_rows = 2 * (max_q + 1) * solution_radius + (max_q + 1) * dim * entry_epsilon
    norm_rows = 2 * (max_q + 1) * dim * solution_radius
    return max(eig_rows, norm_rows)


def certify_node(args, s: mp.mpf, max_q: int, entry_epsilon: mp.mpf):
    mats, vals, e_derivs, lam_derivs, cache_hit = build_center(args, s, max_q)
    dim = mats[0].rows
    vectors = [
        [value / mp.factorial(q) for value in e_derivs[q]]
        for q in range(max_q + 1)
    ]
    lambdas = [lam_derivs[q] / mp.factorial(q) for q in range(max_q + 1)]
    jac = build_jacobian(mats, vectors, lambdas)
    jac_inv = jac ** -1
    inv_norm = mat_inf_norm(jac_inv)
    residual = residual_inf_bound(vectors, entry_epsilon)
    newton_radius = inv_norm * residual
    solution_radius = 2 * newton_radius
    derivative_variation = derivative_variation_bound(
        solution_radius,
        entry_epsilon,
        dim,
        max_q,
    )
    contraction = inv_norm * derivative_variation
    derivative_entry_radius = mp.factorial(max_q) * solution_radius
    derivative_row_l2_radius = mp.sqrt(dim) * derivative_entry_radius
    gap = vals[1] - vals[0] if len(vals) > 1 else mp.inf
    return {
        "s": f(s),
        "cacheHit": cache_hit,
        "gap": gap,
        "jacobianInverseInfNorm": inv_norm,
        "residualInfBound": residual,
        "newtonRadius": newton_radius,
        "solutionRadius": solution_radius,
        "derivativeVariationInfBound": derivative_variation,
        "krawczykContractionBound": contraction,
        "derivativeEntryRadius": derivative_entry_radius,
        "derivativeRowL2Radius": derivative_row_l2_radius,
        "lambda0": vals[0],
        "lambda1": vals[1],
        "maxVectorSeriesCoeffInf": max(vec_inf(row) for row in vectors),
        "maxLambdaSeriesCoeffAbs": max(abs(value) for value in lambdas),
    }


def floats(row: dict) -> dict:
    out = {}
    for key, value in row.items():
        if isinstance(value, mp.mpf):
            out[key] = f(value)
        else:
            out[key] = value
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment-json", default="endpoint_confluent_segment_bernstein_certificate.json")
    parser.add_argument("--trace-tail-json", default="endpoint_confluent_trace_tail_certificate.json")
    parser.add_argument("--kind", choices=["kb"], default="kb")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--condition-matrix-order", type=int, default=90)
    parser.add_argument("--cache-dir", default=".endpoint_eigenrow_interval_cache")
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--needed-trace-q", type=int, default=8)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--node-limit", type=int, default=0)
    parser.add_argument("--trace-radius-target", default="")
    parser.add_argument("--json-out", default="endpoint_eigenrow_interval_propagation.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    segment = load_json(args.segment_json)
    trace_tail = load_json(args.trace_tail_json)
    max_q = args.needed_trace_q
    entry_epsilon = mp.mpf(str(segment["centerTaylorEntryErrorBound"]))
    spectral_epsilon = mp.mpf(str(segment["centerTaylorSpectralErrorBound"]))
    trace_target = (
        mp.mpf(args.trace_radius_target)
        if args.trace_radius_target
        else mp.mpf(str(trace_tail["traceRadiusTarget"]))
    )
    trace_refinement_radius = mp.mpf(str(trace_tail["traceDerivativeEntryRadius"]))
    nodes = collocation_nodes(args)
    if args.node_limit:
        nodes = nodes[: args.node_limit]

    print(
        "Endpoint eigenrow interval propagation "
        f"nodes={len(nodes)} q<={max_q} condition_order={args.condition_matrix_order}",
        flush=True,
    )
    rows = []
    worst = None
    for idx, s in enumerate(nodes):
        row = certify_node(args, s, max_q, entry_epsilon)
        rows.append(row)
        if worst is None or row["derivativeEntryRadius"] > worst["derivativeEntryRadius"]:
            worst = row
        print(
            f"  node {idx+1:02d}/{len(nodes)} s={fmt(s, 8)} "
            f"gap={fmt(row['gap'], 8)} "
            f"K={mp.nstr(row['jacobianInverseInfNorm'], 6)} "
            f"deriv_entry={mp.nstr(row['derivativeEntryRadius'], 6)} "
            f"contr={mp.nstr(row['krawczykContractionBound'], 4)} "
            f"cache={row['cacheHit']}",
            flush=True,
        )

    max_derivative_entry = max(row["derivativeEntryRadius"] for row in rows)
    max_derivative_row = max(row["derivativeRowL2Radius"] for row in rows)
    max_contraction = max(row["krawczykContractionBound"] for row in rows)
    max_newton_radius = max(row["newtonRadius"] for row in rows)
    min_gap = min(row["gap"] for row in rows)
    recurrence_closed = max_contraction < mp.mpf("0.5")
    target_closed = max_derivative_entry < trace_target
    refinement_closed = max_derivative_entry < trace_refinement_radius
    synchronized = args.condition_matrix_order == int(segment["quadratureOrder"])

    data = {
        "theoremName": "endpoint eigenrow interval propagation",
        "segmentJson": args.segment_json,
        "traceTailJson": args.trace_tail_json,
        "conditionMatrixOrder": args.condition_matrix_order,
        "controlledQuadratureOrder": segment["quadratureOrder"],
        "centerSynchronizationPolicy": (
            "The Krawczyk recurrence is conditioned at conditionMatrixOrder. "
            "It is fully synchronized with the deterministic quadrature center "
            "only when conditionMatrixOrder equals controlledQuadratureOrder."
        ),
        "collocationNodeCount": len(nodes),
        "neededTraceQ": max_q,
        "jetOrder": args.jet_order,
        "entryMatrixCoefficientRadius": f(entry_epsilon),
        "spectralMatrixCoefficientRadius": f(spectral_epsilon),
        "traceRadiusTarget": f(trace_target),
        "traceDerivativeRefinementRadius": f(trace_refinement_radius),
        "minGap": f(min_gap),
        "maxNewtonRadius": f(max_newton_radius),
        "maxKrawczykContractionBound": f(max_contraction),
        "maxDerivativeEntryRadius": f(max_derivative_entry),
        "maxDerivativeEntryRadiusText": mp.nstr(max_derivative_entry, 16),
        "maxDerivativeRowL2Radius": f(max_derivative_row),
        "maxDerivativeRowL2RadiusText": mp.nstr(max_derivative_row, 16),
        "derivativeBelowTraceTarget": bool(target_closed),
        "derivativeBelowTraceRefinementRadius": bool(refinement_closed),
        "centerSynchronizedWithSegmentQuadrature": bool(synchronized),
        "worstNode": floats(worst),
        "nodeRows": [floats(row) for row in rows],
        "eigenrowTaylorKrawczykStatus": status(
            "eigenrow Taylor Krawczyk propagation",
            bool(recurrence_closed),
            (
                "Closed when the inverse-Jacobian radius and the local "
                "derivative-variation bound give a contraction for the full "
                "finite eigenvalue/eigenvector Taylor system."
            ),
        ),
        "eigenrowTraceTargetStatus": status(
            "eigenrow derivative interval below trace target",
            bool(target_closed),
            (
                "Closed when the propagated interval for all e^(q), q<=8, is "
                "below the working trace-radius target."
            ),
        ),
        "eigenrowRefinementRadiusStatus": status(
            "eigenrow interval below old refinement radius",
            bool(refinement_closed),
            (
                "Diagnostic only.  The deterministic recurrence bound need "
                "not beat the much smaller empirical 90-point refinement "
                "radius; it must beat the trace target needed by the proof "
                "budget."
            ),
        ),
        "centerSynchronizationStatus": status(
            "recurrence center synchronized with segment quadrature center",
            bool(synchronized),
            (
                "Closed only when the recurrence Jacobian is built at the "
                "same Gauss-Legendre order controlled by the deterministic "
                "Bernstein segment certificate."
            ),
        ),
        "remainingItem": (
            "If centerSynchronizationStatus is open, rerun this certificate "
            "with --condition-matrix-order equal to the controlled quadrature "
            "order, or regenerate the downstream endpoint coefficient "
            "certificate at the displayed synchronized center."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Endpoint eigenrow interval propagation")
    print(f"  max derivative entry radius: {mp.nstr(max_derivative_entry, 12)}")
    print(f"  trace target: {mp.nstr(trace_target, 12)}")
    print(f"  max Krawczyk contraction: {mp.nstr(max_contraction, 12)}")
    print(f"  recurrence closed: {recurrence_closed}")
    print(f"  below trace target: {target_closed}")
    print(f"  center synchronized: {synchronized}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
