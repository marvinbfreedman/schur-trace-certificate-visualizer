#!/usr/bin/env python3
r"""Probe the repair-free Schur complement in trace coordinates.

The current repair-free target is

    Gamma^* Gamma <= C

on the completed trace range.  In a finite Galerkin model this is the
positivity of the Schur complement

    S = C - B^* A^+ B,

where A is the form on ker R, B is the cross block, and C is the trace-side
block.  This script lifts S from the abstract U-basis to the sampled trace
coordinates y = R u.  That lifted matrix is the finite kernel that a continuum
proof should explain.

The output is diagnostic only.  It does not certify the continuum theorem.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp

from quotient_factorization_mp import (
    columns,
    diag,
    frob_norm,
    gram_matrix,
    minmax_eigs,
    positive_part_inverse,
    sampled_centers,
    trace_matrix,
)


def parse_ints(text: str) -> list[int]:
    return [int(piece) for piece in text.replace(",", " ").split()]


def as_float(value) -> float:
    return float(value)


def fmt(value, digits: int = 12) -> str:
    return mp.nstr(value, digits)


def build_args(base: argparse.Namespace, *, basis: int, constraints: int):
    return SimpleNamespace(
        model=base.model,
        kind=base.kind,
        omega=base.omega,
        L=base.L,
        basis=basis,
        quad=max(base.min_quad, base.quad_factor * basis),
        laguerre=base.laguerre,
        endpoint_kernel_order=base.endpoint_kernel_order,
        endpoint_kernel_rmax=base.endpoint_kernel_rmax,
        constraints=constraints,
        constraint_min=base.constraint_min,
        constraint_max=base.constraint_max,
        jet_order=base.jet_order,
        endpoint_order=base.endpoint_order,
        endpoint_rmax=base.endpoint_rmax,
        endpoint_tol=base.endpoint_tol,
        rank_tol=base.rank_tol,
        psd_tol=base.psd_tol,
        margin=base.margin,
        zero_tol=base.zero_tol,
        principal_minor_order=base.principal_minor_order,
        dps=base.dps,
    )


def matrix_to_float(mat: mp.matrix, digits: int | None = None) -> list[list[float]]:
    out = []
    for i in range(mat.rows):
        row = []
        for j in range(mat.cols):
            value = mat[i, j]
            if digits is not None:
                value = mp.mpf(mp.nstr(value, digits))
            row.append(as_float(value))
        out.append(row)
    return out


def normalized_correlation(mat: mp.matrix) -> mp.matrix:
    out = mp.matrix(mat.rows, mat.cols)
    diag_values = [mp.sqrt(max(mp.mpf("0"), mat[i, i])) for i in range(mat.rows)]
    for i in range(mat.rows):
        for j in range(mat.cols):
            denom = diag_values[i] * diag_values[j]
            out[i, j] = mat[i, j] / denom if denom else mp.mpf("0")
    return out


def principal_minor_stats(mat: mp.matrix, max_order: int) -> dict:
    min_det = None
    min_order = None
    negative_count = 0
    checked = 0
    n = mat.rows
    for order in range(1, min(max_order, n) + 1):
        for idx in itertools.combinations(range(n), order):
            sub = mp.matrix(order)
            for i, row in enumerate(idx):
                for j, col in enumerate(idx):
                    sub[i, j] = mat[row, col]
            det = mp.det(sub)
            checked += 1
            if min_det is None or det < min_det:
                min_det = det
                min_order = order
            if det < 0:
                negative_count += 1
    return {
        "checked": checked,
        "maxOrder": min(max_order, n),
        "min": as_float(min_det if min_det is not None else mp.mpf("0")),
        "minOrder": min_order,
        "negativeCount": negative_count,
    }


def contiguous_minor_stats(mat: mp.matrix) -> dict:
    min_det = None
    min_order = None
    negative_count = 0
    checked = 0
    n = mat.rows
    for order in range(1, n + 1):
        for start in range(0, n - order + 1):
            idx = list(range(start, start + order))
            sub = mp.matrix(order)
            for i, row in enumerate(idx):
                for j, col in enumerate(idx):
                    sub[i, j] = mat[row, col]
            det = mp.det(sub)
            checked += 1
            if min_det is None or det < min_det:
                min_det = det
                min_order = order
            if det < 0:
                negative_count += 1
    return {
        "checked": checked,
        "min": as_float(min_det if min_det is not None else mp.mpf("0")),
        "minOrder": min_order,
        "negativeCount": negative_count,
    }


def inverse_structure_stats(mat: mp.matrix, eig_min: mp.mpf) -> dict:
    if mat.rows == 0 or eig_min <= 0:
        return {
            "available": False,
            "reason": "matrix is not strictly positive in the finite model",
        }
    inv = mat**-1
    max_abs = max(abs(inv[i, j]) for i in range(inv.rows) for j in range(inv.cols))
    by_distance = []
    for dist in range(inv.rows):
        values = [
            abs(inv[i, j])
            for i in range(inv.rows)
            for j in range(inv.cols)
            if abs(i - j) == dist
        ]
        if not values:
            continue
        by_distance.append(
            {
                "distance": dist,
                "maxAbs": as_float(max(values)),
                "meanAbs": as_float(mp.fsum(values) / len(values)),
                "maxAbsRelative": as_float(max(values) / max_abs if max_abs else 0),
            }
        )

    diag_nonpositive = sum(1 for i in range(inv.rows) if inv[i, i] <= 0)
    offdiag_positive = sum(
        1
        for i in range(inv.rows)
        for j in range(inv.cols)
        if i != j and inv[i, j] > 0
    )
    return {
        "available": True,
        "frobenius": as_float(frob_norm(inv)),
        "maxAbs": as_float(max_abs),
        "diagNonpositiveCount": diag_nonpositive,
        "offdiagPositiveCount": offdiag_positive,
        "absByIndexDistance": by_distance,
    }


def exponential_mode_fit(centers: list[mp.mpf], values: list[float]) -> dict:
    if not centers or not values or any(value == 0 for value in values):
        return {"available": False}
    x0 = centers[0]
    xs = [center - x0 for center in centers]
    ys = [mp.log(abs(mp.mpf(value))) for value in values]
    n = mp.mpf(len(xs))
    sx = mp.fsum(xs)
    sy = mp.fsum(ys)
    sxx = mp.fsum(x * x for x in xs)
    sxy = mp.fsum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if denom == 0:
        return {"available": False}
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    decay = -slope
    fitted = [mp.e ** (intercept - decay * x) for x in xs]
    target5 = [mp.e ** (-5 * x) for x in xs]
    abs_values = [abs(mp.mpf(value)) for value in values]
    max_fit_residual = max(abs(value - fit) for value, fit in zip(abs_values, fitted))
    max_fit_relative = max(
        abs(value - fit) / max(abs(value), mp.mpf("1e-80"))
        for value, fit in zip(abs_values, fitted)
    )
    max_target5_residual = max(
        abs(value - target) for value, target in zip(abs_values, target5)
    )
    max_target5_relative = max(
        abs(value - target) / max(abs(value), mp.mpf("1e-80"))
        for value, target in zip(abs_values, target5)
    )
    return {
        "available": True,
        "model": "abs(v(a)) ~= exp(intercept-decay*(a-a0))",
        "a0": as_float(x0),
        "decay": as_float(decay),
        "intercept": as_float(intercept),
        "maxFitResidual": as_float(max_fit_residual),
        "maxFitRelativeResidual": as_float(max_fit_relative),
        "maxResidualAgainstExp5": as_float(max_target5_residual),
        "maxRelativeResidualAgainstExp5": as_float(max_target5_relative),
    }


def quotient_blocks(K: mp.matrix, R: mp.matrix, args) -> dict:
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)

    A = N.T * K * N
    B = N.T * K * U
    C = U.T * K * U
    avals, avecs, a_keep, a_zero, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    gamma2 = B.T * Aplus * B if B.cols else mp.matrix(0)
    S = C - gamma2 if B.cols else C

    sigmas = [mp.sqrt(rvals[i]) for i in u_idx]
    RU = R * U
    W = mp.matrix(R.rows, U.cols)
    for j, sigma in enumerate(sigmas):
        for i in range(R.rows):
            W[i, j] = RU[i, j] / sigma
    sigma_inv = diag([1 / sigma for sigma in sigmas]) if sigmas else mp.matrix(0)
    trace_kernel = W * sigma_inv * S * sigma_inv * W.T if sigmas else mp.matrix(R.rows)

    return {
        "rankTol": rank_tol,
        "rvals": rvals,
        "N": N,
        "U": U,
        "A": A,
        "B": B,
        "C": C,
        "Aplus": Aplus,
        "gamma2": gamma2,
        "S": (S + S.T) / 2,
        "traceKernel": (trace_kernel + trace_kernel.T) / 2,
        "sigmas": sigmas,
        "aZeroCount": len(a_zero),
        "aKeepCount": len(a_keep),
    }


def min_mode_trace(blocks: dict, R: mp.matrix) -> list[float]:
    S = blocks["S"]
    if S.rows == 0:
        return []
    vals, vecs = mp.eigsy(S, eigvals_only=False)
    z = columns(vecs, [0])
    y = R * blocks["U"] * z
    scale = max(abs(y[i, 0]) for i in range(y.rows)) if y.rows else mp.mpf("0")
    if scale:
        sign = mp.mpf("1")
        for i in range(y.rows):
            if abs(y[i, 0]) > mp.mpf("1e-80"):
                sign = mp.sign(y[i, 0])
                break
        y = y / (scale * sign)
    return [as_float(y[i, 0]) for i in range(y.rows)]


def probe_one(local) -> dict:
    omega = mp.mpf(local.omega)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, omega, length)
    centers, R = trace_matrix(polys, local)
    blocks = quotient_blocks(K, R, local)

    S = blocks["S"]
    trace_kernel = blocks["traceKernel"]
    corr = normalized_correlation(trace_kernel)

    svals, smin, smax = minmax_eigs(S) if S.rows else ([], mp.mpf("0"), mp.mpf("0"))
    kvals, kmin, kmax = minmax_eigs(K)
    cvals, cmin, cmax = minmax_eigs(blocks["C"]) if blocks["C"].rows else ([], mp.mpf("0"), mp.mpf("0"))
    gvals, gmin, gmax = (
        minmax_eigs(blocks["gamma2"]) if blocks["gamma2"].rows else ([], mp.mpf("0"), mp.mpf("0"))
    )
    tvals, tmin, tmax = (
        minmax_eigs(trace_kernel) if trace_kernel.rows else ([], mp.mpf("0"), mp.mpf("0"))
    )

    diag_values = [trace_kernel[i, i] for i in range(trace_kernel.rows)]
    offdiag_values = [
        trace_kernel[i, j]
        for i in range(trace_kernel.rows)
        for j in range(trace_kernel.cols)
        if i != j
    ]
    principal = principal_minor_stats(trace_kernel, local.principal_minor_order)
    contiguous = contiguous_minor_stats(trace_kernel)
    inverse_stats = inverse_structure_stats(trace_kernel, tmin)
    min_trace = min_mode_trace(blocks, R)
    mode_fit = exponential_mode_fit(centers, min_trace)

    row = {
        "model": local.model,
        "kind": local.kind,
        "omega": as_float(omega),
        "L": as_float(length),
        "basis": local.basis,
        "quad": local.quad,
        "constraints": local.constraints,
        "traceRank": len(blocks["sigmas"]),
        "traceNullity": blocks["N"].cols,
        "aZeroCount": blocks["aZeroCount"],
        "kMin": as_float(kmin),
        "kMax": as_float(kmax),
        "cMin": as_float(cmin),
        "cMax": as_float(cmax),
        "gamma2Min": as_float(gmin),
        "gamma2Max": as_float(gmax),
        "schurComplementMin": as_float(smin),
        "schurComplementMax": as_float(smax),
        "schurComplementEigenvalues": [as_float(v) for v in svals],
        "traceKernelMin": as_float(tmin),
        "traceKernelMax": as_float(tmax),
        "traceKernelEigenvalues": [as_float(v) for v in tvals],
        "traceKernelDiagMin": as_float(min(diag_values) if diag_values else mp.mpf("0")),
        "traceKernelDiagMax": as_float(max(diag_values) if diag_values else mp.mpf("0")),
        "traceKernelOffdiagMin": as_float(min(offdiag_values) if offdiag_values else mp.mpf("0")),
        "traceKernelOffdiagMax": as_float(max(offdiag_values) if offdiag_values else mp.mpf("0")),
        "principalMinorStats": principal,
        "contiguousMinorStats": contiguous,
        "inverseStructureStats": inverse_stats,
        "centers": [as_float(x) for x in centers],
        "minModeTrace": min_trace,
        "minModeExponentialFit": mode_fit,
        "traceKernelMatrix": matrix_to_float(trace_kernel, digits=12),
        "traceKernelCorrelation": matrix_to_float(corr, digits=12),
        "repairFreeFiniteClosed": smin >= -mp.mpf(local.zero_tol),
    }
    row["meaning"] = (
        "S=C-Gamma^*Gamma is the finite repair-free Schur complement. "
        "The traceKernelMatrix is the same form transported to sampled trace "
        "coordinates y=R u."
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="full", choices=["full", "endpoint_b"])
    parser.add_argument("--kind", default="tilde3", choices=["raw1", "raw2", "raw3", "tilde3"])
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis-values", default="6 8 10")
    parser.add_argument("--constraint-values", default="3 5 7")
    parser.add_argument("--quad-factor", type=int, default=2)
    parser.add_argument("--min-quad", type=int, default=14)
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=18)
    parser.add_argument("--endpoint-kernel-rmax", default="10")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=28)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-26")
    parser.add_argument("--margin", default="0")
    parser.add_argument("--zero-tol", default="1e-20")
    parser.add_argument("--principal-minor-order", type=int, default=4)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json-out", default="repair_free_schur_kernel_probe.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    basis_values = parse_ints(args.basis_values)
    constraint_values = parse_ints(args.constraint_values)

    rows = []
    print("repair-free Schur kernel probe", flush=True)
    print(
        "  basis constraints rank K_min S_min S_max trace_kernel_min finite_closed",
        flush=True,
    )
    for basis, constraints in zip(basis_values, constraint_values):
        if constraints >= basis:
            raise SystemExit("each constraints value must be smaller than its basis value")
        local = build_args(args, basis=basis, constraints=constraints)
        row = probe_one(local)
        rows.append(row)
        print(
            f"  {basis:5d} {constraints:11d} {row['traceRank']:4d} "
            f"{fmt(row['kMin'], 8):>12} "
            f"{fmt(row['schurComplementMin'], 8):>12} "
            f"{fmt(row['schurComplementMax'], 8):>12} "
            f"{fmt(row['traceKernelMin'], 8):>16} "
            f"{row['repairFreeFiniteClosed']}",
            flush=True,
        )

    weakest = min(rows, key=lambda row: row["schurComplementMin"]) if rows else None
    data = {
        "theoremName": "repair-free Schur complement trace-kernel probe",
        "meaning": (
            "Finite diagnostic for the continuum target C-Gamma^*Gamma>=0. "
            "Each row transports the Schur complement to sampled trace "
            "coordinates, exposing the finite kernel that a Volterra/Green "
            "proof should factor."
        ),
        "parameters": {
            "model": args.model,
            "kind": args.kind,
            "omega": float(mp.mpf(args.omega)),
            "L": float(mp.mpf(args.L)),
            "basisValues": basis_values,
            "constraintValues": constraint_values,
            "dps": args.dps,
        },
        "rows": rows,
        "weakestRow": weakest,
        "allFiniteRepairFreeClosed": bool(rows) and all(
            row["repairFreeFiniteClosed"] for row in rows
        ),
        "interpretation": [
            "S_min>0 is the finite repair-free condition.",
            (
                "The transported trace kernel is the candidate continuum "
                "operator on X_R; its inverse structure can suggest a "
                "differential/Green proof."
            ),
            (
                "This remains finite evidence until the trace-kernel formula "
                "is derived and bounded in the completed trace range."
            ),
        ],
        "nextProofTarget": (
            "Derive an explicit formula for the trace-side Schur kernel "
            "D_trace=W Sigma^{-1}(C-Gamma^*Gamma)Sigma^{-1}W^*, then prove it "
            "is positive as a Volterra/Green Gram kernel on X_R."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    if weakest:
        print(
            f"  weakest S_min={fmt(weakest['schurComplementMin'], 12)} "
            f"at basis={weakest['basis']} constraints={weakest['constraints']}",
            flush=True,
        )
    print(f"  wrote {args.json_out}", flush=True)


if __name__ == "__main__":
    main()
