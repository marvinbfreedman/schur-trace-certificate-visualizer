#!/usr/bin/env python3
r"""Euler-Lagrange minimizer for the trace-side Schur energy.

For the completed trace map R:V->X_R and N=ker R, choose a section J:X_R->V.
Write f=n+Jx with n in N.  In block form

    Q(n+Jx) = <A n,n> + 2 Re <B x,n> + <C x,x>.

The Euler-Lagrange equation for

    inf { Q(f) : Rf=x }

is

    <A n_x + Bx,h> = 0    for all h in N.

Thus the canonical quotient minimizer is

    n_x = -A^+ Bx,        f_x = Jx - A^+Bx,

with the Moore-Penrose inverse understood in the closed A-form sense.  The
minimized energy is

    Q(f_x) = <(C-B^*A^+B)x,x> = <D_trace x,x>.

This script records the exact Euler-Lagrange theorem and checks the finite
Galerkin version.  It also constructs the finite square representation of the
minimized energy from the positive Schur complement spectrum.

The genuinely open analytic step is stronger: identify the finite spectral
square with an explicit continuum Volterra/Green integral of squares.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp

from quotient_factorization_mp import (
    columns,
    frob_norm,
    gram_matrix,
    positive_part_inverse,
    trace_matrix,
)


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def build_args(args: argparse.Namespace):
    return SimpleNamespace(
        model="full",
        kind="tilde3",
        omega=args.omega,
        L=args.L,
        basis=args.basis,
        quad=max(args.min_quad, args.quad_factor * args.basis),
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=args.constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
        dps=args.dps,
    )


def quotient_blocks(K: mp.matrix, R: mp.matrix, args):
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
    _, _, _, _, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    S = (C - B.T * Aplus * B + (C - B.T * Aplus * B).T) / 2
    return N, U, A, Aplus, B, C, S


def deterministic_vector(size: int, seed: int, scale: str = "13") -> mp.matrix:
    out = mp.matrix(size, 1)
    for i in range(size):
        value = mp.mpf((seed + 2) * (i + 1)) / mp.mpf(scale)
        if (seed + i) % 2:
            value *= -1
        out[i, 0] = value
    return out


def spectral_square(S: mp.matrix, z: mp.matrix, tol: mp.mpf) -> tuple[mp.mpf, mp.mpf, list[float]]:
    vals, vecs = mp.eigsy(S, eigvals_only=False)
    features = []
    total = mp.mpf("0")
    negative_floor = mp.mpf("0")
    coords = vecs.T * z
    for i, val in enumerate(vals):
        if val < -tol:
            negative_floor = min(negative_floor, val)
        lam = max(mp.mpf("0"), val)
        feature = mp.sqrt(lam) * coords[i, 0]
        features.append(float(feature))
        total += feature * feature
    return total, negative_floor, features


def finite_minimizer_check(args: argparse.Namespace) -> dict:
    local = build_args(args)
    omega = mp.mpf(local.omega)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, omega, length)
    _, R = trace_matrix(polys, local)
    N, U, A, Aplus, B, C, S = quotient_blocks(K, R, local)

    svals = mp.eigsy(S, eigvals_only=True)
    rows = []
    max_el_residual = mp.mpf("0")
    max_trace_error = mp.mpf("0")
    max_energy_error = mp.mpf("0")
    max_square_error = mp.mpf("0")
    min_perturbation_gap = None
    negative_square_floor = mp.mpf("0")

    for seed in range(args.test_vectors):
        z = deterministic_vector(U.cols, seed)
        n = -Aplus * B * z
        f = N * n + U * z
        trace_target = R * U * z
        trace_error = frob_norm(R * f - trace_target)
        el_residual = frob_norm(A * n + B * z)
        energy = (f.T * K * f)[0]
        schur_energy = (z.T * S * z)[0]
        square_energy, neg_floor, features = spectral_square(S, z, mp.mpf(args.psd_tol))
        negative_square_floor = min(negative_square_floor, neg_floor)

        max_trace_error = max(max_trace_error, trace_error)
        max_el_residual = max(max_el_residual, el_residual)
        max_energy_error = max(max_energy_error, abs(energy - schur_energy))
        max_square_error = max(max_square_error, abs(square_energy - schur_energy))

        perturbation_gaps = []
        for pseed in range(3):
            h = deterministic_vector(N.cols, seed + pseed + 7, scale="29")
            fp = f + N * h
            gap = (fp.T * K * fp)[0] - energy
            predicted = (h.T * A * h)[0]
            perturbation_gaps.append(float(gap))
            if min_perturbation_gap is None or gap < min_perturbation_gap:
                min_perturbation_gap = gap
            max_energy_error = max(max_energy_error, abs(gap - predicted))

        rows.append(
            {
                "seed": seed,
                "traceError": float(trace_error),
                "eulerLagrangeResidual": float(el_residual),
                "energy": float(energy),
                "schurEnergy": float(schur_energy),
                "spectralSquareEnergy": float(square_energy),
                "energyMinusSchur": float(energy - schur_energy),
                "squareMinusSchur": float(square_energy - schur_energy),
                "perturbationGaps": perturbation_gaps,
                "firstFeatures": features[: min(6, len(features))],
            }
        )

    closed = (
        max_el_residual < mp.mpf(args.check_tol)
        and max_trace_error < mp.mpf(args.check_tol)
        and max_energy_error < mp.mpf(args.check_tol)
        and max_square_error < mp.mpf(args.check_tol)
        and (min_perturbation_gap is None or min_perturbation_gap >= -mp.mpf(args.check_tol))
        and negative_square_floor >= -mp.mpf(args.check_tol)
    )

    return {
        "basis": local.basis,
        "constraints": local.constraints,
        "traceRank": U.cols,
        "traceNullity": N.cols,
        "schurComplementMin": float(svals[0]) if len(svals) else 0.0,
        "schurComplementMax": float(svals[-1]) if len(svals) else 0.0,
        "maxEulerLagrangeResidual": float(max_el_residual),
        "maxTraceError": float(max_trace_error),
        "maxEnergyIdentityError": float(max_energy_error),
        "maxSpectralSquareError": float(max_square_error),
        "minPerturbationGap": float(min_perturbation_gap or mp.mpf("0")),
        "negativeSpectralFloor": float(negative_square_floor),
        "rows": rows,
        "closed": closed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--basis", type=int, default=8)
    parser.add_argument("--constraints", type=int, default=5)
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
    parser.add_argument("--check-tol", default="1e-35")
    parser.add_argument("--test-vectors", type=int, default=4)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="trace_euler_lagrange_minimizer.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    finite = finite_minimizer_check(args)

    statuses = {
        "eulerLagrangeEquationStatus": status(
            "Euler-Lagrange equation on trace fiber",
            True,
            (
                "For fixed trace x, variations h in ker R give "
                "Q(f_x,h)=0, equivalently A n_x + Bx=0."
            ),
        ),
        "minimizerFormulaStatus": status(
            "Moore-Penrose minimizer formula",
            True,
            (
                "The canonical minimizer is f_x=Jx-A^+Bx.  Other minimizers "
                "can differ only by A-null vectors in ker R."
            ),
        ),
        "minimizedEnergyStatus": status(
            "minimized energy is Schur complement",
            True,
            "Substitution gives Q(f_x)=<x,(C-B^*A^+B)x>.",
        ),
        "finiteMinimizerCheckStatus": status(
            "finite Euler-Lagrange minimizer check",
            finite["closed"],
            (
                "Finite Galerkin minimizers satisfy the trace constraint, "
                "Euler-Lagrange residual, energy identity, and perturbation "
                "minimality checks below tolerance."
            ),
            blocker=None if finite["closed"] else "A finite minimizer identity failed.",
        ),
        "finiteSpectralSquareStatus": status(
            "finite minimized energy square representation",
            finite["closed"],
            (
                "Since the finite Schur complement is positive, its spectral "
                "decomposition writes the minimized energy as a finite sum of "
                "squares."
            ),
            blocker=None if finite["closed"] else "The finite square identity failed.",
        ),
        "continuumVolterraGreenSquareStatus": status(
            "continuum Volterra/Green integral square",
            False,
            (
                "The minimizer and abstract square are solved, but the square "
                "has not yet been identified with explicit Volterra/Green "
                "features for the continuum kernel."
            ),
            blocker=(
                "Derive the Green operator solving A n=-Bx in the continuum "
                "and prove C-B^*A^+B equals an integral of explicit residual "
                "Green features."
            ),
        ),
    }

    data = {
        "theoremName": "trace-fiber Euler-Lagrange minimizer",
        "exactEquations": [
            "f=n+Jx, n in N=ker R",
            "Q(n+Jx)=<An,n>+2 Re <Bx,n>+<Cx,x>",
            "Euler-Lagrange: A n_x + Bx = 0 in N^*",
            "canonical minimizer: f_x=Jx-A^+Bx",
            "minimized energy: Q(f_x)=<x,(C-B^*A^+B)x>",
            (
                "finite spectral square: <Sz,z>=sum_k "
                "|sqrt(lambda_k)<v_k,z>|^2 when S>=0"
            ),
        ],
        "greenInterpretation": {
            "operatorEquation": "A n_x = -B x",
            "affineTraceConstraint": "R(Jx+n_x)=x",
            "variationalBoundaryCondition": "Q(f_x,h)=0 for every h with Rh=0",
            "dirichletToNeumannForm": "D_trace(x,y)=Q(f_x,f_y)",
            "squareTarget": (
                "Find explicit residual Green features G_x(u) such that "
                "D_trace(x,y)=integral G_x(u)G_y(u)dmu(u)."
            ),
        },
        "statuses": statuses,
        "finiteCheck": finite,
        "continuumTheoremClosed": False,
        "correctedConclusion": (
            "The Euler-Lagrange minimizer is solved exactly and verified in "
            "finite quotient coordinates.  The minimized energy is the Schur "
            "kernel and has a finite spectral square in every positive finite "
            "probe.  The open part is the explicit continuum Volterra/Green "
            "integral-square identification."
        ),
        "nextProofTarget": (
            "Build the continuum Green solver for A n=-Bx from the Volterra "
            "Sturm operator, then compute the residual feature map whose "
            "L2/Volterra norm equals C-B^*A^+B."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("trace-fiber Euler-Lagrange minimizer")
    print(f"  finite minimizer check: {finite['closed']}")
    print(f"  max EL residual: {finite['maxEulerLagrangeResidual']:.3e}")
    print(f"  max trace error: {finite['maxTraceError']:.3e}")
    print(f"  max energy identity error: {finite['maxEnergyIdentityError']:.3e}")
    print(f"  finite Schur min: {finite['schurComplementMin']:.12e}")
    print("  continuum Volterra/Green square: open")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
