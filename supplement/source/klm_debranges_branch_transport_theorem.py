#!/usr/bin/env python3
r"""Exact and closed-cone branch transport criteria.

This is the next bridge layer after the canonical Hardy image

    K_E = <h^+,h^+> - <h^-,h^->.

There are two possible ways to finish the KLM-to-de Branges bridge.

Exact route.
    Construct a single isometry/unitary U such that

        U h_z^+ = G_+(z),       U h_z^- = G_-(z).

    This is equivalent to equality of the joint Gram kernels of the family
    {h_z^+, h_z^-} and the proposed Volterra/KLM feature family
    {G_+(z), G_-(z)}.  At present the missing object is not an abstract
    Hilbert-space theorem; it is the concrete map z -> x_z into the completed
    Volterra trace image.

Closed-cone route.
    Truncate the Hardy features:

        h_{z,R}^\pm = 1_{[0,R]} h_z^\pm.

    These converge strongly in L^2(0,infty), and their signed Gram kernels
    converge entrywise to K_E.  Therefore, if each truncated signed Gram is a
    KLM pullback or a positive-cone limit of KLM pullbacks, then K_E is positive
    by the endpoint closed-cone theorem.  The truncation convergence is proved
    here; the remaining missing piece is again the Volterra/KLM approximation
    map for the truncated Hardy branches.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path

from klm_debranges_pullback_probe import XiTransform, eigvalsh_hermitian, frobenius


def default_z_nodes() -> list[complex]:
    return [
        0.0 + 0.35j,
        0.45 + 0.40j,
        0.95 + 0.55j,
        1.55 + 0.75j,
        2.25 + 0.95j,
    ]


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


def hardy_feature_scalar(xi: XiTransform, omega: float, z: complex, branch: str) -> complex:
    if branch == "plus":
        return xi(z + 1j * omega) / math.sqrt(2.0 * math.pi)
    if branch == "minus":
        return xi(z - 1j * omega) / math.sqrt(2.0 * math.pi)
    raise ValueError(branch)


def halfline_inner(a_z: complex, a_w: complex, z: complex, w: complex, rmax: float | None) -> complex:
    delta = z - w.conjugate()
    if rmax is None:
        integral = 1.0 / (1j * (w.conjugate() - z))
    else:
        integral = (cmath.exp(1j * delta * rmax) - 1.0) / (1j * delta)
    return a_z * a_w.conjugate() * integral


def joint_hardy_gram(xi: XiTransform, omega: float, nodes: list[complex], rmax: float | None):
    labels = []
    for branch in ("plus", "minus"):
        for idx, z in enumerate(nodes):
            labels.append((branch, idx, z, hardy_feature_scalar(xi, omega, z, branch)))
    n = len(labels)
    gram = [[0j for _ in range(n)] for _ in range(n)]
    for i, (_bi, _ii, w, aw) in enumerate(labels):
        for j, (_bj, _jj, z, az) in enumerate(labels):
            gram[i][j] = halfline_inner(az, aw, z, w, rmax)
    return gram


def signed_kernel_from_joint(joint: list[list[complex]], count: int) -> list[list[complex]]:
    return [[joint[i][j] - joint[i + count][j + count] for j in range(count)] for i in range(count)]


def matrix_sub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a))] for i in range(len(a))]


def hardy_tail_bound(xi: XiTransform, omega: float, nodes: list[complex], rmax: float) -> float:
    max_bound = 0.0
    for w in nodes:
        for z in nodes:
            denom = z.imag + w.imag
            plus = abs(xi(z + 1j * omega) * xi(w + 1j * omega).conjugate())
            minus = abs(xi(z - 1j * omega) * xi(w - 1j * omega).conjugate())
            bound = (plus + minus) * math.exp(-denom * rmax) / (2.0 * math.pi * denom)
            max_bound = max(max_bound, bound)
    return max_bound


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xi-tmax", type=float, default=8.0)
    parser.add_argument("--xi-intervals", type=int, default=1400)
    parser.add_argument("--rmax", type=float, default=40.0)
    parser.add_argument("--json-out", default="klm_debranges_branch_transport_theorem.json")
    args = parser.parse_args()

    xi = XiTransform(args.xi_tmax, args.xi_intervals)
    nodes = default_z_nodes()
    count = len(nodes)
    joint_full = joint_hardy_gram(xi, args.omega, nodes, None)
    joint_trunc = joint_hardy_gram(xi, args.omega, nodes, args.rmax)
    signed_full = signed_kernel_from_joint(joint_full, count)
    signed_trunc = signed_kernel_from_joint(joint_trunc, count)
    signed_error = matrix_sub(signed_full, signed_trunc)
    joint_error = matrix_sub(joint_full, joint_trunc)
    signed_rel = frobenius(signed_error) / max(frobenius(signed_full), 1e-300)
    joint_rel = frobenius(joint_error) / max(frobenius(joint_full), 1e-300)
    entry_tail_bound = hardy_tail_bound(xi, args.omega, nodes, args.rmax)

    exact_transport_criterion = {
        "statement": (
            "A single isometry U with U h_z^+=G_+(z) and U h_z^-=G_-(z) exists "
            "on the finite evaluation span iff the joint Gram matrix of "
            "{h_z^+,h_z^-} equals the joint Gram matrix of {G_+(z),G_-(z)}. "
            "The continuum statement is the corresponding Kolmogorov/isometry "
            "extension after density."
        ),
        "requiredJointGramEqualities": [
            "<h_z^+,h_w^+> = <G_+(z),G_+(w)>",
            "<h_z^-,h_w^-> = <G_-(z),G_-(w)>",
            "<h_z^+,h_w^-> = <G_+(z),G_-(w)>",
            "<h_z^-,h_w^+> = <G_-(z),G_+(w)>",
        ],
        "missingConcreteMap": "z -> x_z in the completed Volterra trace image, with G_±(z)=G_±(f_{x_z})",
    }

    data = {
        "theoremName": "exact and strong-limit branch transport theorem",
        "omega": args.omega,
        "xiQuadrature": {"tmax": args.xi_tmax, "intervals": args.xi_intervals},
        "rmax": args.rmax,
        "zNodes": [{"real": z.real, "imag": z.imag} for z in nodes],
        "statuses": {
            "exactUCriterionStatus": status(
                "exact branch transport criterion",
                True,
                (
                    "The exact U problem is reduced to equality of the joint "
                    "Hardy and Volterra/KLM feature Gram kernels.  This is the "
                    "standard Hilbert-space isometry extension criterion."
                ),
            ),
            "exactUConstructedStatus": status(
                "exact branch transport U constructed",
                False,
                (
                    "The criterion is known, but no concrete z -> x_z map into "
                    "the completed Volterra trace image has been constructed."
                ),
                blocker=exact_transport_criterion["missingConcreteMap"],
            ),
            "hardyStrongTruncationStatus": status(
                "Hardy branch strong truncation",
                True,
                (
                    "h_{z,R}^±=1_[0,R]h_z^± converges strongly to h_z^± in "
                    "L^2(0,infty).  The finite signed Gram error is explicitly "
                    "bounded by the exponential tail displayed in this JSON."
                ),
            ),
            "klmClosedConeLimitStatus": status(
                "KLM closed-cone branch limit",
                False,
                (
                    "The Hardy truncation limit is closed, but the truncated "
                    "Hardy branches have not yet been represented as KLM "
                    "pullbacks or strong limits of KLM/Volterra feature branches."
                ),
                blocker=(
                    "Construct x_{z,R} or phase-space pullback vectors such that "
                    "G_±(x_{z,R}) -> h_{z,R}^± jointly, with Gram convergence."
                ),
            ),
        },
        "exactTransportCriterion": exact_transport_criterion,
        "hardyStrongLimitProof": {
            "truncatedFeature": "h_{z,R}^± = 1_[0,R](r) (2*pi)^(-1/2) E_omega^±(z) exp(i z r)",
            "tailNormSquared": "|E_omega^±(z)|^2 exp(-2 Im(z) R)/(4*pi*Im(z))",
            "entryTailBoundMax": entry_tail_bound,
            "jointGramRelativeTailOnSample": joint_rel,
            "signedKernelRelativeTailOnSample": signed_rel,
            "signedKernelTailEigenMin": min(eigvalsh_hermitian(signed_error)),
            "signedKernelTailEigenMax": max(eigvalsh_hermitian(signed_error)),
        },
        "closedConeBridgeClosed": False,
        "nextProofTarget": (
            "Construct the concrete evaluation trace map z -> x_z, or truncated "
            "maps x_{z,R}, so that the Volterra/KLM branch features have the "
            "same joint Gram limits as the Hardy branches."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Exact and strong-limit branch transport theorem")
    print("  exact U criterion: closed")
    print("  exact U constructed: False")
    print("  Hardy strong truncation: closed")
    print(f"  sample joint Gram truncation rel error: {joint_rel:.6e}")
    print(f"  sample signed kernel truncation rel error: {signed_rel:.6e}")
    print(f"  max entry tail bound: {entry_tail_bound:.6e}")
    print("  KLM closed-cone branch limit: False")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
