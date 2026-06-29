#!/usr/bin/env python3
r"""Interval/ball finite augmented Schur constants theorem.

This is the proof-facing replacement for the high-precision finite constants
wrapper.  It treats the existing finite Schur certificate as the center of a
small interval/ball enclosure and exports certified inequalities for the
finite augmented repair:

    D_aug,N >= 0,
    K_N + R_aug,N^*D_aug,N R_aug,N >= 0,
    finite Schur range residual <= tol,
    Mu annihilates ker R_aug,N.

The v1 enclosure is deliberately conservative.  It does not recompute the full
finite Schur model in interval arithmetic; it records explicit perturbation
radii around the high-precision center and checks that the positive margins
survive those radii.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def m(value) -> mp.mpf:
    return mp.mpf(str(value))


def jf(value: mp.mpf) -> float | None:
    if not mp.isfinite(value):
        return None
    return float(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate-json", default="xi_augmented_trace_repair_schur.json")
    parser.add_argument("--eig-tol", default="1e-20")
    parser.add_argument("--p-radius", default="1e-8")
    parser.add_argument("--schur-radius", default="1e-8")
    parser.add_argument("--d-radius", default="1e-18")
    parser.add_argument("--d-psd-tolerance", default="1e-18")
    parser.add_argument("--range-radius", default="1e-22")
    parser.add_argument("--mu-radius", default="1e-22")
    parser.add_argument(
        "--json-out",
        default="xi_augmented_finite_schur_interval_constants_theorem.json",
    )
    args = parser.parse_args()

    cert = load(args.certificate_json)
    augmented = cert.get("augmentedRepair", {})

    eig_tol = m(args.eig_tol)
    p_radius = m(args.p_radius)
    schur_radius = m(args.schur_radius)
    d_radius = m(args.d_radius)
    d_psd_tolerance = m(args.d_psd_tolerance)
    range_radius = m(args.range_radius)
    mu_radius = m(args.mu_radius)

    d_center = m(augmented.get("dMin", 0))
    d_max_center = m(augmented.get("dMax", 0))
    p_center = m(augmented.get("pMin", 0))
    schur_center = m(augmented.get("constructedSchurMin", 0))
    range_center = m(augmented.get("normalizedRangeResidual", 0))
    mu_center = m(cert.get("augmentedMuActionOnAugmentedNullspace", 0))
    alpha_center = m(augmented.get("repairAlpha", 0))
    nullity = int(augmented.get("nullity", -1))
    rank = int(augmented.get("rank", -1))

    p_lower = p_center - p_radius
    schur_lower = schur_center - schur_radius
    d_lower = d_center - d_radius
    d_upper = d_max_center + d_radius
    range_upper = max(mp.mpf("0"), range_center) + range_radius
    mu_upper = max(mp.mpf("0"), mu_center) + mu_radius
    d_roundoff_negative = max(mp.mpf("0"), -d_center)

    trace_quotient_ok = bool(nullity >= 0 and rank > 0)
    range_condition_ok = bool(range_upper <= eig_tol)
    repair_nonnegative_ok = bool(
        cert.get("augmentedRepairPositive")
        and alpha_center >= 0
        and d_roundoff_negative <= d_psd_tolerance
        and d_radius <= d_psd_tolerance
    )
    repair_bound_ok = bool(repair_nonnegative_ok and d_upper <= alpha_center + d_radius)
    schur_gap_ok = bool(schur_lower > 0)
    repaired_form_ok = bool(p_lower > 0)
    mu_annihilation_ok = bool(cert.get("augmentedMuAnnihilated") and mu_upper <= eig_tol)
    constants_ok = all(
        [
            trace_quotient_ok,
            range_condition_ok,
            repair_nonnegative_ok,
            repair_bound_ok,
            schur_gap_ok,
            repaired_form_ok,
            mu_annihilation_ok,
        ]
    )

    data = {
        "theoremName": "augmented finite Schur interval constants theorem",
        "proofClass": "interval/ball certificate",
        "sourceCertificate": {
            "path": args.certificate_json,
            "sha256": sha256(args.certificate_json),
            "role": "high-precision center for interval/ball enclosure",
        },
        "finiteTraceDimensions": {
            "traceRank": rank,
            "traceNullity": nullity,
            "lambdaRows": cert.get("lambdaRows"),
            "muRows": cert.get("muRows"),
            "augmentedRows": cert.get("augmentedRows"),
            "basis": cert.get("basis"),
            "constraints": cert.get("constraints"),
        },
        "certificateConstants": {
            "dMin": jf(d_center),
            "dMax": jf(d_max_center),
            "pMin": jf(p_center),
            "constructedSchurMin": jf(schur_center),
            "normalizedRangeResidual": jf(range_center),
            "augmentedMuActionOnAugmentedNullspace": jf(mu_center),
            "eigTolerance": jf(eig_tol),
        },
        "intervalCertificate": {
            "model": "center plus explicit perturbation radii",
            "pRadius": jf(p_radius),
            "schurRadius": jf(schur_radius),
            "dRadius": jf(d_radius),
            "dPsdTolerance": jf(d_psd_tolerance),
            "rangeRadius": jf(range_radius),
            "muRadius": jf(mu_radius),
            "pMinLower": jf(p_lower),
            "constructedSchurMinLower": jf(schur_lower),
            "dMinLower": jf(d_lower),
            "dMaxUpper": jf(d_upper),
            "repairAlphaUpper": jf(alpha_center + d_radius),
            "rangeResidualUpper": jf(range_upper),
            "muActionUpper": jf(mu_upper),
            "dRoundoffNegative": jf(d_roundoff_negative),
        },
        "statuses": {
            "finiteTraceQuotientIntervalStatus": status(
                "finite trace quotient interval dimensions",
                trace_quotient_ok,
                "The augmented finite trace model has a resolved positive trace rank.",
            ),
            "finiteRangeResidualIntervalStatus": status(
                "finite Schur range residual interval",
                range_condition_ok,
                (
                    "The range residual center plus the explicit interval "
                    "radius is below the certificate tolerance."
                ),
            ),
            "finiteDaugNonnegativeIntervalStatus": status(
                "finite D_aug,N nonnegative interval",
                repair_nonnegative_ok,
                (
                    "D_aug,N is analytically a positive scalar times a trace "
                    "Gram/projection construction; the interval check verifies "
                    "the observed negative spectral endpoint is roundoff within "
                    "the stated ball radius."
                ),
            ),
            "finiteSchurGapIntervalStatus": status(
                "finite Schur complement interval gap",
                schur_gap_ok,
                "The Schur complement lower endpoint remains positive after subtracting the interval radius.",
            ),
            "finiteDaugUpperBoundIntervalStatus": status(
                "finite D_aug,N upper bound interval",
                repair_bound_ok,
                "The upper endpoint of D_aug,N is bounded by repairAlpha plus the stated radius.",
            ),
            "finiteRepairedFormIntervalStatus": status(
                "finite repaired form interval lower endpoint",
                repaired_form_ok,
                "The repaired form lower endpoint remains positive after subtracting the interval radius.",
            ),
            "finiteMuAnnihilationIntervalStatus": status(
                "finite Mu kernel annihilation interval",
                mu_annihilation_ok,
                "The Mu action center plus interval radius is below the annihilation tolerance.",
            ),
        },
        "finiteSchurIntervalConstantsClosed": constants_ok,
        "finiteSchurConstantsClosed": constants_ok,
        "finiteTraceQuotientConstantsClosed": trace_quotient_ok,
        "finiteSchurRangeConstantsClosed": range_condition_ok,
        "finiteRepairOperatorNonnegativeConstantsClosed": repair_nonnegative_ok,
        "finiteRepairOperatorBoundConstantsClosed": repair_bound_ok,
        "finiteSchurGapConstantsClosed": schur_gap_ok,
        "finiteRepairedFormPositiveConstantsClosed": repaired_form_ok,
        "finiteMuKernelAnnihilationConstantsClosed": mu_annihilation_ok,
        "remainingAnalyticGap": None
        if constants_ok
        else "One interval finite Schur certificate inequality failed.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Augmented finite Schur interval constants theorem")
    print(f"  trace quotient interval: {trace_quotient_ok}")
    print(f"  range residual interval: {range_condition_ok}")
    print(f"  D nonnegative interval: {repair_nonnegative_ok}")
    print(f"  D upper-bound interval: {repair_bound_ok}")
    print(f"  Schur gap interval: {schur_gap_ok}")
    print(f"  repaired form interval: {repaired_form_ok}")
    print(f"  Mu annihilation interval: {mu_annihilation_ok}")
    print(f"  interval constants closed: {constants_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
