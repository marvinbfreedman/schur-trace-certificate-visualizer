#!/usr/bin/env python3
r"""Summarize the repair-free route after primitive boundary collapse.

Inputs:
  * primitive_boundary_transport_audit.json proves D_bdy=0;
  * primitive_trace_image_density.json proves the primitive trace image is
    dense in X_R;
  * dq_vanishing_schur_defect_scan*.json give finite tests of
    Gamma^*Gamma-C <= 0.

Conclusion:
  The remaining theorem is the continuum operator inequality

      Gamma^*Gamma-C <= 0 on X_R.

If proved, D_q=0 and the original Weyl lift follows from the exact primitive
transport without any boundary repair.  The present script records finite
evidence only; it does not close the continuum inequality.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_existing(paths: list[str]) -> list[dict]:
    out = []
    for path in paths:
        candidate = Path(path)
        if candidate.exists():
            data = load(path)
            data["_path"] = path
            out.append(data)
    return out


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primitive-boundary-json",
        default="primitive_boundary_transport_audit.json",
    )
    parser.add_argument(
        "--primitive-density-json",
        default="primitive_trace_image_density.json",
    )
    parser.add_argument(
        "--scan-jsons",
        default=(
            "dq_vanishing_schur_defect_scan.json "
            "dq_vanishing_schur_defect_scan_b8c5.json "
            "dq_vanishing_schur_defect_scan_b10c7.json"
        ),
    )
    parser.add_argument(
        "--primitive-endpoint-json",
        default="primitive_endpoint_compatibility_theorem.json",
    )
    parser.add_argument("--json-out", default="dq_vanishing_repair_route_summary.json")
    args = parser.parse_args()

    primitive = load(args.primitive_boundary_json)
    density = load(args.primitive_density_json)
    primitive_endpoint = load(args.primitive_endpoint_json)
    scans = load_existing(args.scan_jsons.replace(",", " ").split())

    primitive_zero = bool(
        primitive.get("statuses", {})
        .get("canonicalPrimitiveBoundaryZeroStatus", {})
        .get("closed")
    )
    trace_dense = bool(density.get("primitiveTraceImageDenseInXR"))
    dq_equiv = bool(density.get("dqZeroOnPrimitiveImageEquivalentToDqZero"))
    endpoint_dq_zero = bool(
        primitive_endpoint.get("statuses", {})
        .get("dqVanishesOnXRStatus", {})
        .get("closed")
    )

    scan_rows = []
    for scan in scans:
        for row in scan.get("rows", []):
            item = dict(row)
            item["sourceJson"] = scan["_path"]
            scan_rows.append(item)

    all_finite_zero = bool(scan_rows) and all(row.get("finiteDqZero") for row in scan_rows)
    worst = max(scan_rows, key=lambda row: float(row["schurDefectMax"])) if scan_rows else None

    primitive_status = status(
        "primitive boundary repair vanishes",
        primitive_zero,
        (
            "The primitive transport has no lower or upper endpoint term, so "
            "D_bdy=0."
        ),
    )
    density_status = status(
        "primitive trace image dense in X_R",
        trace_dense and dq_equiv,
        (
            "The primitive image contains the compact smooth Volterra core; "
            "after completion its trace image is all of X_R.  Thus D_q|Y=0 "
            "is equivalent to D_q=0 on X_R."
        ),
    )
    finite_status = status(
        "finite Schur-defect scans support D_q=0",
        all_finite_zero,
        (
            (
                "All scanned finite models have lambda_max(Gamma^*Gamma-C)<0, "
                "so the finite positive part D_q is zero."
            )
            if all_finite_zero
            else "A scanned finite model has positive Schur defect."
        ),
        blocker=None
        if all_finite_zero
        else "Inspect the positive finite Schur-defect mode.",
    )
    continuum_status = status(
        "continuum repair-free theorem",
        endpoint_dq_zero,
        (
            (
                "Closed by the primitive endpoint compatibility theorem: the "
                "Green-lift contraction gives M<=P on the completed trace "
                "image, so D_q=(M-P)_+=0, equivalently "
                "Gamma^*Gamma-C<=0 on X_R."
            )
            if endpoint_dq_zero
            else (
                "The finite evidence is negative, but the continuum inequality "
                "Gamma^*Gamma-C<=0 on X_R has not been proved."
            )
        ),
        blocker=None
        if endpoint_dq_zero
        else (
            "Prove the continuum anti-Douglas inequality "
            "Gamma^*Gamma<=C on X_R, or derive an equivalent Gram/Volterra "
            "factorization showing Q_Phi>=0 directly."
        ),
    )

    data = {
        "theoremName": "D_q-vanishing repair route summary",
        "primitiveBoundaryJson": args.primitive_boundary_json,
        "primitiveDensityJson": args.primitive_density_json,
        "primitiveEndpointJson": args.primitive_endpoint_json,
        "scanJsons": [scan["_path"] for scan in scans],
        "statuses": {
            "primitiveBoundaryZeroStatus": primitive_status,
            "primitiveTraceDensityStatus": density_status,
            "finiteSchurDefectStatus": finite_status,
            "continuumRepairFreeTheoremStatus": continuum_status,
        },
        "scanRows": scan_rows,
        "worstFiniteRow": worst,
        "allFiniteDqZero": all_finite_zero,
        "continuumTheoremClosed": continuum_status["closed"],
        "correctedNextProofTarget": continuum_status["blocker"],
        "formalImplicationIfClosed": (
            "If Gamma^*Gamma-C<=0 on X_R, then D_q=0.  Since D_bdy=0 and the "
            "primitive transport is exact, Q_original=Q_Phi=||G_q F||^2 plus "
            "a nonnegative Schur block, giving original Weyl positivity in "
            "the normalized full-Phi model."
        ),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("D_q vanishing repair route summary")
    print(f"  primitive D_bdy=0: {primitive_zero}")
    print(f"  primitive trace dense: {trace_dense}")
    print(f"  D_q|Y=0 iff D_q=0: {dq_equiv}")
    print(f"  scan rows: {len(scan_rows)}")
    if worst:
        print(
            "  worst finite h_max="
            f"{worst['schurDefectMax']:.12e} "
            f"basis={worst['basis']} constraints={worst['constraints']}"
        )
    print(f"  all finite D_q zero: {all_finite_zero}")
    print(f"  continuum theorem closed: {continuum_status['closed']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
