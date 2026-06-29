#!/usr/bin/env python3
r"""Ball certificate for full active row rank of the endpoint map.

The exact endpoint Green BVP has active endpoint equation

    M z = -b(d).

If the active endpoint map M has full row rank, then ``ker M^*={0}`` and
Fredholm compatibility is vacuous.  This script certifies that conclusion from
an explicit ball enclosure for M:

    M_exact in {M0 + E : ||E||_F <= delta}.

Then Weyl's singular-value perturbation bound gives

    sigma_min(M_exact) >= sigma_min(M0) - delta.

The script is intentionally strict about logical status.  A positive lower
bound proves full row rank only relative to the supplied entry ball.  It is a
fully rigorous continuum certificate only when that ball has been produced by
a rigorous interval/ODE enclosure; otherwise it is a conditional or diagnostic
rank certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


def load_json(path: str) -> dict:
    candidate = Path(path)
    if not candidate.exists():
        raise SystemExit(f"missing {path}")
    return json.loads(candidate.read_text(encoding="utf-8"))


def f(x) -> float:
    if not mp.isfinite(x):
        return None
    return float(x)


def fmt(x, digits: int = 12) -> str:
    return mp.nstr(x, digits)


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def mp_matrix(rows: list[list[float | int | str]]) -> mp.matrix:
    return mp.matrix([[mp.mpf(str(value)) for value in row] for row in rows])


def frobenius(mat: mp.matrix) -> mp.mpf:
    return mp.sqrt(
        mp.fsum(abs(mat[i, j]) ** 2 for i in range(mat.rows) for j in range(mat.cols))
    )


def vector_norm(values) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in values))


def gram_eigenvalues_from_json(data: dict, mat: mp.matrix) -> list[mp.mpf]:
    stored = data.get("endpointMapGramEigenvalues")
    if stored:
        return sorted(mp.mpf(str(value)) for value in stored)
    gram = (mat * mat.T + mat * mat.T) / 2
    return sorted(mp.eigsy(gram, eigvals_only=True))


def entry_radii(center: mp.matrix, rel_radius: mp.mpf, abs_radius: mp.mpf) -> mp.matrix:
    radii = mp.matrix(center.rows, center.cols)
    for i in range(center.rows):
        for j in range(center.cols):
            radii[i, j] = abs_radius + rel_radius * abs(center[i, j])
    return radii


def add_variant_radii(radii: mp.matrix, center: mp.matrix, variant_paths: list[str]) -> list[dict]:
    summaries = []
    for path in variant_paths:
        data = load_json(path)
        if "endpointMap" not in data:
            raise SystemExit(f"{path} has no endpointMap")
        variant = mp_matrix(data["endpointMap"])
        if variant.rows != center.rows or variant.cols != center.cols:
            raise SystemExit(f"{path} endpointMap shape mismatch")
        max_delta = mp.mpf("0")
        frob_delta2 = mp.mpf("0")
        for i in range(center.rows):
            for j in range(center.cols):
                delta = abs(variant[i, j] - center[i, j])
                radii[i, j] = max(radii[i, j], delta)
                max_delta = max(max_delta, delta)
                frob_delta2 += delta * delta
        summaries.append(
            {
                "path": path,
                "maxEntryDelta": f(max_delta),
                "frobeniusDelta": f(mp.sqrt(frob_delta2)),
            }
        )
    return summaries


def max_entry(mat: mp.matrix) -> mp.mpf:
    return max(
        [abs(mat[i, j]) for i in range(mat.rows) for j in range(mat.cols)]
        or [mp.mpf("0")]
    )


def least_left_direction(data: dict) -> list[mp.mpf]:
    stored = data.get("leftObstructionVector")
    if stored:
        return [mp.mpf(str(value)) for value in stored]
    mat = mp_matrix(data["endpointMap"])
    gram = (mat * mat.T + mat * mat.T) / 2
    vals, vecs = mp.eigsy(gram, eigvals_only=False)
    idx = min(range(len(vals)), key=lambda i: abs(vals[i]))
    vec = [vecs[i, idx] for i in range(mat.rows)]
    nrm = vector_norm(vec)
    if nrm:
        vec = [value / nrm for value in vec]
    if vec and vec[0] < 0:
        vec = [-value for value in vec]
    return vec


def source_pairings(data: dict, direction: list[mp.mpf]) -> dict:
    rows = []
    max_abs = mp.mpf("0")
    max_rel = mp.mpf("0")
    for row in data.get("rows", []):
        endpoint_vector = row.get("endpointVector")
        if not endpoint_vector:
            continue
        values = [mp.mpf(str(value)) for value in endpoint_vector]
        pairing = mp.fsum(direction[i] * values[i] for i in range(len(values)))
        norm = vector_norm(values)
        rel = abs(pairing) / max(mp.mpf("1"), norm)
        max_abs = max(max_abs, abs(pairing))
        max_rel = max(max_rel, rel)
        rows.append(
            {
                "sourceNode": row.get("sourceNode"),
                "component": row.get("component"),
                "endpointVectorNorm": f(norm),
                "leftDirectionPairing": f(pairing),
                "leftDirectionRelativePairing": f(rel),
            }
        )
    return {
        "maxAbsPairing": f(max_abs),
        "maxRelativePairing": f(max_rel),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-json", default="adjoint_green_endpoint_selection.json")
    parser.add_argument(
        "--variant-json",
        action="append",
        default=[],
        help="Additional endpoint JSON files used to enlarge the empirical entry ball.",
    )
    parser.add_argument("--entry-rel-radius", default="1e-22")
    parser.add_argument("--entry-abs-radius", default="0")
    parser.add_argument("--frobenius-radius", default="0")
    parser.add_argument(
        "--assume-rigorous-entry-ball",
        action="store_true",
        help="Mark the supplied entry/frobenius ball as a rigorous enclosure of the exact endpoint map.",
    )
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--json-out", default="endpoint_map_rank_ball_certificate.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    data = load_json(args.endpoint_json)
    if "endpointMap" not in data:
        raise SystemExit(
            f"{args.endpoint_json} has no endpointMap; rerun adjoint_green_endpoint_selection.py"
        )

    center = mp_matrix(data["endpointMap"])
    eigs = gram_eigenvalues_from_json(data, center)
    positive_eigs = [max(mp.mpf("0"), value) for value in eigs]
    sigma_min_center = mp.sqrt(min(positive_eigs))
    sigma_max_center = mp.sqrt(max(positive_eigs))
    center_frob = mp.mpf(str(data.get("endpointMapFrobenius", frobenius(center))))
    rel_radius = mp.mpf(args.entry_rel_radius)
    abs_radius = mp.mpf(args.entry_abs_radius)
    radii = entry_radii(center, rel_radius, abs_radius)
    variant_summaries = add_variant_radii(radii, center, args.variant_json)
    entry_ball_frob = frobenius(radii)
    explicit_frob = mp.mpf(args.frobenius_radius)
    total_delta = max(entry_ball_frob, explicit_frob)
    lower_sigma = max(mp.mpf("0"), sigma_min_center - total_delta)
    lower_gram = lower_sigma * lower_sigma
    full_rank_conditional = bool(lower_sigma > 0)
    rigorous_full_rank = bool(full_rank_conditional and args.assume_rigorous_entry_ball)

    active_dim = int(data.get("activeDim", center.rows))
    cols = center.cols
    required_frob_radius = sigma_min_center
    required_uniform_entry_radius = sigma_min_center / mp.sqrt(active_dim * cols)
    required_relative_frob_radius = (
        sigma_min_center / center_frob if center_frob else mp.inf
    )
    condition_upper = (
        (sigma_max_center + total_delta) / lower_sigma if lower_sigma else mp.inf
    )
    pinv_upper = 1 / lower_sigma if lower_sigma else mp.inf
    left_dir = least_left_direction(data)
    pairings = source_pairings(data, left_dir)

    conditional_status = status(
        "conditional ball full active row rank",
        full_rank_conditional,
        (
            "Weyl perturbation gives sigma_min(M_exact)>=sigma_min(M0)-delta "
            "for the supplied endpoint-map ball."
        ),
    )
    rigorous_status = status(
        "rigorous endpoint map full active row rank",
        rigorous_full_rank,
        (
            "Closed only if the supplied endpoint-map ball is a rigorous "
            "interval enclosure of the exact fundamental-matrix endpoint map."
        ),
    )
    kernel_status = status(
        "ker M^* is zero",
        rigorous_full_rank,
        (
            "If the rigorous full-row-rank status is closed, then the left "
            "Fredholm obstruction space is trivial and endpoint compatibility "
            "is vacuous."
        ),
    )
    pinv_status = status(
        "Moore-Penrose endpoint bound",
        full_rank_conditional,
        (
            "On the supplied ball, ||M^+|| <= 1/(sigma_min(M0)-delta).  This "
            "becomes a rigorous continuum bound only with a rigorous endpoint "
            "map enclosure."
        ),
    )

    data_out = {
        "theoremName": "endpoint map rank ball certificate",
        "endpointJson": args.endpoint_json,
        "variantJson": args.variant_json,
        "radiusModel": {
            "entryRelRadius": args.entry_rel_radius,
            "entryAbsRadius": args.entry_abs_radius,
            "frobeniusRadius": args.frobenius_radius,
            "assumeRigorousEntryBall": args.assume_rigorous_entry_ball,
            "interpretation": (
                "The rank conclusion is rigorous only if these radii enclose "
                "the exact endpoint fundamental-matrix map."
            ),
        },
        "activeDim": active_dim,
        "endpointColumns": cols,
        "centerGramEigenvalues": [f(value) for value in eigs],
        "centerSigmaMin": f(sigma_min_center),
        "centerSigmaMax": f(sigma_max_center),
        "centerConditionNumber": f(sigma_max_center / sigma_min_center),
        "centerFrobeniusNorm": f(center_frob),
        "centerMaxEntry": f(max_entry(center)),
        "entryBallFrobeniusRadius": f(entry_ball_frob),
        "totalFrobeniusRadius": f(total_delta),
        "lowerSigmaMin": f(lower_sigma),
        "lowerGramEigenvalue": f(lower_gram),
        "pinvNormUpperBound": f(pinv_upper),
        "conditionNumberUpperBound": f(condition_upper),
        "maxAllowedFrobeniusRadiusForRank": f(required_frob_radius),
        "maxAllowedUniformEntryRadiusForRank": f(required_uniform_entry_radius),
        "maxAllowedRelativeFrobeniusRadiusForRank": f(required_relative_frob_radius),
        "variantSummaries": variant_summaries,
        "leastStableLeftDirection": [f(value) for value in left_dir],
        "leastStableSourcePairings": pairings,
        "conditionalBallFullRowRankStatus": conditional_status,
        "rigorousEndpointMapFullRowRankStatus": rigorous_status,
        "kernelMStarZeroStatus": kernel_status,
        "moorePenroseEndpointBoundStatus": pinv_status,
        "proof": [
            "Let M_exact=M0+E with ||E||_2<=||E||_F<=delta.",
            "Weyl's singular-value inequality gives sigma_min(M_exact)>=sigma_min(M0)-delta.",
            "If this lower bound is positive, M_exact has full active row rank.",
            "Full active row rank implies ker(M^*)={0}, so the endpoint Fredholm compatibility condition is vacuous.",
            "The same lower bound gives ||M^+||<=1/(sigma_min(M0)-delta).",
        ],
        "remainingAnalyticGap": (
            "Produce a rigorous interval enclosure for the exact endpoint "
            "fundamental-matrix map with Frobenius radius below the displayed "
            "maxAllowedFrobeniusRadiusForRank."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data_out, indent=2) + "\n", encoding="utf-8")

    print("Endpoint map rank ball certificate")
    print(f"  center sigma_min: {fmt(sigma_min_center, 12)}")
    print(f"  center sigma_max: {fmt(sigma_max_center, 12)}")
    print(f"  relative singular margin: {fmt(required_relative_frob_radius, 12)}")
    print(f"  supplied Frobenius radius: {fmt(total_delta, 12)}")
    print(f"  lower sigma_min: {fmt(lower_sigma, 12)}")
    print(f"  conditional full row rank: {full_rank_conditional}")
    print(f"  rigorous full row rank: {rigorous_full_rank}")
    print(f"  ||M^+|| upper: {fmt(pinv_upper, 12)}")
    print(f"  left-direction max relative source pairing: {pairings['maxRelativePairing']:.12e}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
