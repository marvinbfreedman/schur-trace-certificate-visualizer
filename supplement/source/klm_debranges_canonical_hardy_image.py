#!/usr/bin/env python3
r"""Canonical Hardy image for the shifted-Xi de Branges kernel.

For

    E_omega(z)  = Xi(z+i omega),
    E_omega#(z) = Xi(z-i omega),

the de Branges kernel has the exact half-line Hardy feature identity

    K_E(w,z)
      = <h_z^+, h_w^+>_{L^2(0,infty)}
        - <h_z^-, h_w^->_{L^2(0,infty)},

with

    h_z^+(r) = (2 pi)^(-1/2) E_omega(z)  exp(i z r),
    h_z^-(r) = (2 pi)^(-1/2) E_omega#(z) exp(i z r).

Indeed, for Im z, Im w > 0,

    integral_0^infty exp(i(z-conj(w))r) dr = 1/(i(conj(w)-z)).

This is the canonical de Branges evaluation-vector image in the Fourier/Hardy
model.  It is not yet the missing KLM pullback.  The remaining bridge is to
identify these two Hardy branches with the completed Weyl/Volterra/KLM
positive and negative branches, or to prove that the signed Hardy form is a
closed positive-cone limit of KLM pullbacks.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path

from klm_debranges_pullback_probe import XiTransform, debranges_kernel, eigvalsh_hermitian, frobenius


def default_z_nodes() -> list[complex]:
    return [
        0.0 + 0.35j,
        0.45 + 0.40j,
        0.95 + 0.55j,
        1.55 + 0.75j,
        2.25 + 0.95j,
    ]


def hardy_branch_inner(
    xi: XiTransform,
    omega: float,
    w: complex,
    z: complex,
    *,
    branch: str,
) -> complex:
    if branch == "plus":
        ez = xi(z + 1j * omega)
        ew = xi(w + 1j * omega)
    elif branch == "minus":
        ez = xi(z - 1j * omega)
        ew = xi(w - 1j * omega)
    else:
        raise ValueError(branch)
    return ez * ew.conjugate() / (2j * math.pi * (w.conjugate() - z))


def hardy_branch_inner_truncated(
    xi: XiTransform,
    omega: float,
    w: complex,
    z: complex,
    *,
    branch: str,
    rmax: float,
) -> complex:
    if branch == "plus":
        ez = xi(z + 1j * omega)
        ew = xi(w + 1j * omega)
    elif branch == "minus":
        ez = xi(z - 1j * omega)
        ew = xi(w - 1j * omega)
    else:
        raise ValueError(branch)
    delta = z - w.conjugate()
    integral = (cmath.exp(1j * delta * rmax) - 1.0) / (1j * delta)
    return ez * ew.conjugate() * integral / (2.0 * math.pi)


def matrix_stats(mat: list[list[complex]]) -> dict:
    eigs = eigvalsh_hermitian(mat)
    return {
        "frobenius": frobenius(mat),
        "eigenMin": min(eigs),
        "eigenMax": max(eigs),
        "eigenvalues": eigs,
    }


def generalized_branch_eigs(minus: list[list[complex]], plus: list[list[complex]]) -> dict:
    try:
        import numpy as np

        p = np.array(plus, dtype=complex)
        m = np.array(minus, dtype=complex)
        p = 0.5 * (p + p.conjugate().T)
        m = 0.5 * (m + m.conjugate().T)
        vals, vecs = np.linalg.eigh(p)
        cutoff = max(float(vals[-1]) * 1e-12, 1e-18)
        keep = vals > cutoff
        if not bool(np.any(keep)):
            return {"available": False, "reason": "plus branch has no resolved positive spectrum"}
        inv_sqrt = (vecs[:, keep] / np.sqrt(vals[keep])) @ vecs[:, keep].conjugate().T
        compressed = inv_sqrt @ m @ inv_sqrt
        compressed = 0.5 * (compressed + compressed.conjugate().T)
        ceigs = [float(x) for x in np.linalg.eigvalsh(compressed)]
        rank = int(np.count_nonzero(keep))
    except Exception:
        compressed = cholesky_congruence(minus, plus)
        ceigs = eigvalsh_hermitian(compressed)
        cutoff = 0.0
        rank = len(plus)
    return {
        "available": True,
        "plusEigenCutoff": cutoff,
        "resolvedRank": rank,
        "maxGeneralizedEigenvalue": max(ceigs),
        "minGeneralizedEigenvalue": min(ceigs),
        "eigenvalues": ceigs,
        "sampleContractionHolds": max(ceigs) <= 1.0 + 1e-10,
        "sampleMargin": 1.0 - max(ceigs),
    }


def cholesky_congruence(minus: list[list[complex]], plus: list[list[complex]]) -> list[list[complex]]:
    """Return L^{-1} minus L^{-*} for plus=L L^*."""
    n = len(plus)
    lower = [[0j for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = plus[i][j]
            for k in range(j):
                s -= lower[i][k] * lower[j][k].conjugate()
            if i == j:
                lower[i][j] = math.sqrt(max(s.real, 0.0))
            else:
                lower[i][j] = s / lower[j][j]

    def solve_lower(rhs: list[complex]) -> list[complex]:
        out = [0j for _ in range(n)]
        for i in range(n):
            s = rhs[i]
            for k in range(i):
                s -= lower[i][k] * out[k]
            out[i] = s / lower[i][i]
        return out

    # A = L^{-1} minus L^{-*}; compute columns by applying L^{-1} to minus
    # times columns of L^{-*}.
    inv_l = []
    for j in range(n):
        e = [0j for _ in range(n)]
        e[j] = 1.0
        inv_l.append(solve_lower(e))

    # inv_l_cols[j][i] is (L^{-1})_{ij}; build dense inverse.
    linv = [[inv_l[j][i] for j in range(n)] for i in range(n)]
    tmp = [[0j for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            tmp[i][j] = sum(linv[i][k] * minus[k][j] for k in range(n))
    out = [[0j for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            out[i][j] = sum(tmp[i][k] * linv[j][k].conjugate() for k in range(n))
    return out


def subtract(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    return [[a[i][j] - b[i][j] for j in range(len(a))] for i in range(len(a))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=1400)
    parser.add_argument("--rmax", type=float, default=40.0)
    parser.add_argument("--json-out", default="klm_debranges_canonical_hardy_image.json")
    args = parser.parse_args()

    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    nodes = default_z_nodes()
    n = len(nodes)
    direct = [[0j for _ in range(n)] for _ in range(n)]
    hardy_exact = [[0j for _ in range(n)] for _ in range(n)]
    hardy_truncated = [[0j for _ in range(n)] for _ in range(n)]
    plus = [[0j for _ in range(n)] for _ in range(n)]
    minus = [[0j for _ in range(n)] for _ in range(n)]

    for i, w in enumerate(nodes):
        for j, z in enumerate(nodes):
            direct[i][j] = debranges_kernel(xi, args.omega, w, z)
            plus[i][j] = hardy_branch_inner(xi, args.omega, w, z, branch="plus")
            minus[i][j] = hardy_branch_inner(xi, args.omega, w, z, branch="minus")
            hardy_exact[i][j] = plus[i][j] - minus[i][j]
            hardy_truncated[i][j] = (
                hardy_branch_inner_truncated(xi, args.omega, w, z, branch="plus", rmax=args.rmax)
                - hardy_branch_inner_truncated(xi, args.omega, w, z, branch="minus", rmax=args.rmax)
            )

    exact_error = subtract(direct, hardy_exact)
    trunc_error = subtract(direct, hardy_truncated)
    exact_rel = frobenius(exact_error) / max(frobenius(direct), 1e-300)
    trunc_rel = frobenius(trunc_error) / max(frobenius(direct), 1e-300)

    data = {
        "theoremName": "canonical Hardy image of shifted-Xi de Branges kernel",
        "omega": args.omega,
        "xiQuadrature": {"tmax": args.xi_tmax, "intervals": args.xi_intervals},
        "rmax": args.rmax,
        "zNodes": [{"real": z.real, "imag": z.imag} for z in nodes],
        "formula": {
            "Eplus": "E_omega(z)=Xi(z+i omega)",
            "Eminus": "E_omega#(z)=Xi(z-i omega)",
            "hPlus": "(2*pi)^(-1/2) E_omega(z) exp(i z r), r>=0",
            "hMinus": "(2*pi)^(-1/2) E_omega#(z) exp(i z r), r>=0",
            "kernelIdentity": "K_E(w,z)=<h_z^+,h_w^+>-<h_z^-,h_w^->",
            "halfLineIntegral": "int_0^infty exp(i(z-conj(w))r)dr=1/(i(conj(w)-z))",
        },
        "directKernelStats": matrix_stats(direct),
        "plusBranchStats": matrix_stats(plus),
        "minusBranchStats": matrix_stats(minus),
        "hardyBranchContraction": generalized_branch_eigs(minus, plus),
        "hardyExactErrorFrobenius": frobenius(exact_error),
        "hardyExactRelativeError": exact_rel,
        "hardyTruncatedErrorFrobenius": frobenius(trunc_error),
        "hardyTruncatedRelativeError": trunc_rel,
        "canonicalHardyImageClosed": exact_rel < 1e-12,
        "directKlmPullbackClosed": False,
        "closedPositiveConeLimitStatus": {
            "closed": False,
            "reducedTo": (
                "Identify h_z^+ and h_z^- with the completed Weyl/Volterra/KLM "
                "positive and negative branch features, or prove they are strong "
                "limits of KLM pullback features.  Then the already proved "
                "branch contraction gives K_E >= 0 by closed-cone closure."
            ),
        },
        "nextProofTarget": (
            "Construct a unitary or isometric transport U from the Hardy branch "
            "features h_z^± to the Volterra/KLM branch features G_±, and prove "
            "U h_z^- = C K E U h_z^+ on the closed evaluation span."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Canonical Hardy image for shifted-Xi de Branges kernel")
    print(f"  omega={args.omega:g} z_nodes={n}")
    print(f"  exact Hardy identity relative error: {exact_rel:.6e}")
    print(f"  truncated r<={args.rmax:g} relative error: {trunc_rel:.6e}")
    print(f"  target min eigenvalue: {data['directKernelStats']['eigenMin']:.6e}")
    print(f"  plus branch min eigenvalue: {data['plusBranchStats']['eigenMin']:.6e}")
    print(f"  minus branch min eigenvalue: {data['minusBranchStats']['eigenMin']:.6e}")
    if data["hardyBranchContraction"].get("available"):
        print(
            "  max eig(plus^-1/2 minus plus^-1/2): "
            f"{data['hardyBranchContraction']['maxGeneralizedEigenvalue']:.12e}"
        )
        print(f"  sample contraction margin: {data['hardyBranchContraction']['sampleMargin']:.12e}")
    print(f"  canonical Hardy image closed: {data['canonicalHardyImageClosed']}")
    print("  direct KLM pullback closed: False")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
