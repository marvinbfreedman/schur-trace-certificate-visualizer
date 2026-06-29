#!/usr/bin/env python3
r"""Probe candidate KLM-to-de Branges pullbacks.

The completed KLM theorem gives a positive kernel on phase-space points

    K_KLM(p,q) = Q_omega(p-q) exp(i sigma(p,q)/2),   p=(s,t).

The missing bridge would produce an explicit map T_omega from de Branges
evaluation vectors to the KLM positive cone, ideally

    K_Eomega(w,z) = <T k_w, K_KLM T k_z>.

This script tests finite pullback ansatzes of that form.  It is a diagnostic,
not a theorem: it builds vectors alpha_z on a phase-space grid and compares

    alpha_w^* K_KLM alpha_z

with the shifted-Xi de Branges kernel.  If a candidate succeeds, the residual
matrix should be small or positive semidefinite after one scalar normalization.
If it fails, the residual records the next missing conjugation/weight/limit.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path

from klm_test import SymplecticFourier, phi

try:
    import numpy as np
except ImportError:  # pragma: no cover - the repo has pure-python fallbacks elsewhere.
    np = None


def simpson_grid_symmetric(xmax: float, intervals: int) -> tuple[list[float], list[float]]:
    if intervals % 2:
        intervals += 1
    h = 2.0 * xmax / intervals
    xs = [-xmax + i * h for i in range(intervals + 1)]
    ws = []
    for i in range(intervals + 1):
        if i == 0 or i == intervals:
            c = 1.0
        elif i % 2:
            c = 4.0
        else:
            c = 2.0
        ws.append(c * h / 3.0)
    return xs, ws


class XiTransform:
    def __init__(self, tmax: float, intervals: int):
        self.ts, self.ws = simpson_grid_symmetric(tmax, intervals)

    def __call__(self, z: complex) -> complex:
        total = 0j
        for t, w in zip(self.ts, self.ws):
            total += w * phi(t) * cmath.exp(1j * z * t)
        return total


def debranges_kernel(xi: XiTransform, omega: float, w: complex, z: complex) -> complex:
    ez = xi(z + 1j * omega)
    ew = xi(w + 1j * omega)
    esz = xi(z - 1j * omega)
    esw = xi(w - 1j * omega)
    return (ez * ew.conjugate() - esz * esw.conjugate()) / (2j * math.pi * (w.conjugate() - z))


def phase_space_grid(smax: float, tmax: float, ns: int, nt: int) -> list[tuple[float, float, float]]:
    points = []
    for i in range(ns):
        s = -smax + 2.0 * smax * i / (ns - 1) if ns > 1 else 0.0
        for j in range(nt):
            t = -tmax + 2.0 * tmax * j / (nt - 1) if nt > 1 else 0.0
            points.append((s, t, 1.0))
    return points


def klm_kernel_matrix(points: list[tuple[float, float, float]], qfun: SymplecticFourier):
    n = len(points)
    mat = [[0j for _ in range(n)] for _ in range(n)]
    for i, (si, ti, _) in enumerate(points):
        for j, (sj, tj, _) in enumerate(points):
            q = qfun(si - sj, ti - tj)
            phase = cmath.exp(0.5j * (ti * sj - si * tj))
            mat[i][j] = q * phase
    return mat


def coherent_vector(
    z: complex,
    points: list[tuple[float, float, float]],
    *,
    kind: str,
    width: float,
    omega: float,
) -> list[complex]:
    x = z.real
    y = z.imag
    out = []
    for s, t, weight in points:
        gaussian = math.exp(-0.5 * ((s - x) ** 2 + (t - y) ** 2) / (width * width))
        if kind == "bargmann":
            phase = cmath.exp(0.5j * (t * x - s * y))
            branch = math.exp(-omega * s)
        elif kind == "anti_bargmann":
            phase = cmath.exp(-0.5j * (t * x - s * y))
            branch = math.exp(omega * s)
        elif kind == "weyl_plane":
            phase = cmath.exp(1j * (t * x - s * y))
            branch = 1.0
        elif kind == "shifted_plane":
            phase = cmath.exp(1j * (t * x - s * y))
            branch = math.exp(-omega * s)
        else:
            raise ValueError(kind)
        out.append(math.sqrt(weight) * gaussian * branch * phase)
    return out


def quadratic_pullback(alpha_w: list[complex], kk: list[list[complex]], alpha_z: list[complex]) -> complex:
    total = 0j
    n = len(alpha_z)
    for i in range(n):
        row_total = 0j
        for j in range(n):
            row_total += kk[i][j] * alpha_z[j]
        total += alpha_w[i].conjugate() * row_total
    return total


def eigvalsh_hermitian(mat: list[list[complex]]) -> list[float]:
    if np is not None:
        arr = np.array(mat, dtype=complex)
        arr = 0.5 * (arr + arr.conjugate().T)
        return [float(x) for x in np.linalg.eigvalsh(arr)]
    # Small pure-python fallback via realification.
    from klm_test import hermitian_real_form, jacobi_eigenvalues_symmetric

    return jacobi_eigenvalues_symmetric(hermitian_real_form(mat))


def frobenius(mat: list[list[complex]]) -> float:
    return math.sqrt(sum(abs(v) ** 2 for row in mat for v in row))


def fit_scalar(target: list[list[complex]], pull: list[list[complex]]) -> complex:
    num = 0j
    den = 0.0
    n = len(target)
    for i in range(n):
        for j in range(n):
            num += pull[i][j].conjugate() * target[i][j]
            den += abs(pull[i][j]) ** 2
    return num / den if den else 0j


def matrix_for_candidate(args, z_nodes, points, kk, xi, kind: str):
    alphas = [
        coherent_vector(z, points, kind=kind, width=args.width, omega=args.omega)
        for z in z_nodes
    ]
    n = len(z_nodes)
    pull = [[0j for _ in range(n)] for _ in range(n)]
    target = [[0j for _ in range(n)] for _ in range(n)]
    for i, w in enumerate(z_nodes):
        for j, z in enumerate(z_nodes):
            pull[i][j] = quadratic_pullback(alphas[i], kk, alphas[j])
            target[i][j] = debranges_kernel(xi, args.omega, w, z)
    scale = fit_scalar(target, pull)
    residual = [[target[i][j] - scale * pull[i][j] for j in range(n)] for i in range(n)]
    target_eigs = eigvalsh_hermitian(target)
    pull_eigs = eigvalsh_hermitian([[scale * v for v in row] for row in pull])
    residual_eigs = eigvalsh_hermitian(residual)
    return {
        "kind": kind,
        "scaleReal": scale.real,
        "scaleImag": scale.imag,
        "targetFrobenius": frobenius(target),
        "scaledPullbackFrobenius": frobenius([[scale * v for v in row] for row in pull]),
        "residualFrobenius": frobenius(residual),
        "relativeResidualFrobenius": frobenius(residual) / max(frobenius(target), 1e-300),
        "targetEigenMin": min(target_eigs),
        "targetEigenMax": max(target_eigs),
        "scaledPullbackEigenMin": min(pull_eigs),
        "scaledPullbackEigenMax": max(pull_eigs),
        "residualEigenMin": min(residual_eigs),
        "residualEigenMax": max(residual_eigs),
        "residualLooksPositive": min(residual_eigs) >= -1e-10 * max(1.0, max(abs(x) for x in residual_eigs)),
        "residualLooksZero": frobenius(residual) <= 1e-8 * max(frobenius(target), 1e-300),
    }


def default_z_nodes() -> list[complex]:
    return [
        0.0 + 0.35j,
        0.45 + 0.40j,
        0.95 + 0.55j,
        1.55 + 0.75j,
        2.25 + 0.95j,
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--q-ymax", type=float, default=8.0)
    parser.add_argument("--q-intervals", type=int, default=700)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=1400)
    parser.add_argument("--smax", type=float, default=3.0)
    parser.add_argument("--tmax", type=float, default=3.0)
    parser.add_argument("--ns", type=int, default=9)
    parser.add_argument("--nt", type=int, default=9)
    parser.add_argument("--width", type=float, default=0.85)
    parser.add_argument(
        "--width-scan",
        default="",
        help="Comma-separated Gaussian widths. If set, overrides --width for candidate scanning.",
    )
    parser.add_argument("--json-out", default="klm_debranges_pullback_probe.json")
    args = parser.parse_args()

    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    qfun = SymplecticFourier(args.omega, args.q_ymax, args.q_intervals)
    points = phase_space_grid(args.smax, args.tmax, args.ns, args.nt)
    kk = klm_kernel_matrix(points, qfun)
    z_nodes = default_z_nodes()
    widths = (
        [float(x) for x in args.width_scan.split(",") if x.strip()]
        if args.width_scan
        else [args.width]
    )
    candidates = []
    original_width = args.width
    for width in widths:
        args.width = width
        for kind in ("bargmann", "anti_bargmann", "weyl_plane", "shifted_plane"):
            item = matrix_for_candidate(args, z_nodes, points, kk, xi, kind)
            item["width"] = width
            candidates.append(item)
    args.width = original_width
    best = min(candidates, key=lambda item: item["relativeResidualFrobenius"])
    residual_psd = any(item["residualLooksPositive"] for item in candidates)
    residual_zero = any(item["residualLooksZero"] for item in candidates)

    data = {
        "theoremName": "KLM-to-de Branges finite pullback probe",
        "omega": args.omega,
        "phaseSpaceGrid": {
            "smax": args.smax,
            "tmax": args.tmax,
            "ns": args.ns,
            "nt": args.nt,
            "count": len(points),
        },
        "xiQuadrature": {"tmax": args.xi_tmax, "intervals": args.xi_intervals},
        "qQuadrature": {"ymax": args.q_ymax, "intervals": args.q_intervals},
        "widths": widths,
        "zNodes": [{"real": z.real, "imag": z.imag} for z in z_nodes],
        "candidates": candidates,
        "bestCandidate": best,
        "foundExactFinitePullback": residual_zero,
        "foundPositiveResidualDomination": residual_psd,
        "diagnosis": (
            "This finite coherent-grid ansatz did not find the exact bridge; "
            "the next analytic target is to replace the ad hoc Gaussian vectors "
            "by the Weyl image of the de Branges evaluation functional or to "
            "derive a boundary/weight correction from the residual."
            if not residual_zero
            else "A finite candidate residual is numerically zero; derive the exact formula."
        ),
        "nextProofTarget": (
            "Derive the canonical Weyl/Bargmann image of k_z from the Fourier "
            "representation Xi(z)=int Phi(t) exp(i z t) dt, then rerun this "
            "probe with that non-ad-hoc vector."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("KLM-to-de Branges finite pullback probe")
    print(f"  omega={args.omega:g} phase grid={len(points)} z_nodes={len(z_nodes)}")
    for item in candidates:
        print(
            "  {kind:14s} width={width:.3g} rel_res={rel:.6e} res_min={rmin:.6e} "
            "res_max={rmax:.6e} scale={sr:.6e}{si:+.2e}i".format(
                kind=item["kind"],
                width=item["width"],
                rel=item["relativeResidualFrobenius"],
                rmin=item["residualEigenMin"],
                rmax=item["residualEigenMax"],
                sr=item["scaleReal"],
                si=item["scaleImag"],
            )
        )
    print(f"  best={best['kind']} rel_res={best['relativeResidualFrobenius']:.6e}")
    print(f"  exact finite pullback found: {residual_zero}")
    print(f"  positive residual domination found: {residual_psd}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
