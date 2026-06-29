#!/usr/bin/env python3
"""Build the JSON data used by visual_explainer.html.

The page is meant as an educational dashboard, not a proof certificate.  It
uses the same endpoint B-model jet computations as the research scripts, but
keeps the sample count modest so it can be regenerated quickly.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from endpoint_defect_family_mp import defect_at  # noqa: E402
from endpoint_kb_confluent_mp import integrate  # noqa: E402


def f(value):
    return float(value)


def signed(vec):
    if vec[0] > 0:
        return [-x for x in vec]
    return vec


def matrix_to_float(mat):
    return [[f(mat[i, j]) for j in range(mat.cols)] for i in range(mat.rows)]


def main() -> None:
    mp.mp.dps = 50
    n = 9
    rmax = mp.mpf("12")
    order = 35
    tol = mp.mpf("1e-20")
    s0 = mp.mpf("0.5")

    jet, _ = integrate("kb", mp.pi, s0, n, rmax, order)
    vals, vecs = mp.eigsy(jet, eigvals_only=False)
    lowest = signed([vecs[i, 0] for i in range(n)])

    centers = [mp.mpf("0.02") + i * (mp.mpf("0.545") - mp.mpf("0.02")) / 7 for i in range(8)]
    defect_rows = []
    prev = None
    for center in centers:
        row_vals, neg, vec = defect_at(center, n, rmax, order, tol)
        vec = signed(vec)
        if prev is not None and mp.fsum(a * b for a, b in zip(prev, vec)) < 0:
            vec = [-x for x in vec]
        prev = vec
        defect_rows.append(
            {
                "s": f(center),
                "lambda0": f(row_vals[0]),
                "lambda1": f(row_vals[1]),
                "negCount": len(neg),
                "coeffs": [f(x) for x in vec],
                "pivot": max(range(n), key=lambda k: abs(vec[k])),
            }
        )

    data = {
        "meta": {
            "title": "Endpoint Trace / Schur Visualizer",
            "generatedBy": "build_visual_data.py",
            "dps": 50,
            "jetOrder": n,
            "quadratureOrder": order,
            "rmax": f(rmax),
            "center": f(s0),
            "sStar": 0.5530736,
            "e8Zero": 0.1646167,
        },
        "jet": {
            "matrix": matrix_to_float(jet),
            "eigenvalues": [f(x) for x in vals],
            "lowestEigenvector": [f(x) for x in lowest],
        },
        "defectFamily": defect_rows,
        "quotient": {
            "endpointB": {
                "label": "Endpoint B stress model",
                "domain": "[0,2]",
                "basis": 16,
                "constraints": 8,
                "kMin": -1.7980435e-13,
                "kerMin": 1.8483932e-10,
                "crossNorm": 2.3951005e-5,
                "rangeResidual": 0.0,
                "gamma2Max": 6.9003097e-8,
                "schurDefectMax": 1.8209307e-13,
                "repairedMin": 9.8314105e-19,
            },
            "tilde3": {
                "label": "Full tilde3 reduced kernel",
                "domain": "[0,8]",
                "basis": 12,
                "constraints": 6,
                "kMin": 9.3829204e-7,
                "kerMin": 2.5040653e-5,
                "crossNorm": 9.2643019e-2,
                "rangeResidual": 0.0,
                "gamma2Max": 5.8883244e-2,
                "schurDefectMax": -2.0143861e-6,
                "repairedMin": 9.3829204e-7,
            },
        },
        "douglasScan": [
            {
                "basis": 10,
                "constraints": 5,
                "kMin": 2.3022135e-12,
                "kerMin": 2.6504592e-6,
                "rangeResidual": 0.0,
                "gamma2Max": 1.3580683e-6,
                "repairedMin": 2.3022145e-12,
            },
            {
                "basis": 12,
                "constraints": 6,
                "kMin": -1.7246722e-13,
                "kerMin": 7.9617137e-8,
                "rangeResidual": 0.0,
                "gamma2Max": 5.2114947e-7,
                "repairedMin": 9.9964651e-19,
            },
            {
                "basis": 14,
                "constraints": 7,
                "kMin": -1.7952118e-13,
                "kerMin": 6.8426087e-9,
                "rangeResidual": 0.0,
                "gamma2Max": 3.0884044e-7,
                "repairedMin": 9.9929372e-19,
            },
            {
                "basis": 16,
                "constraints": 8,
                "kMin": -1.7980435e-13,
                "kerMin": 1.8483932e-10,
                "rangeResidual": 0.0,
                "gamma2Max": 6.9003097e-8,
                "repairedMin": 9.8314105e-19,
            },
        ],
        "explain": {
            "jetMatrix": (
                "A jet matrix records all derivatives of a kernel at one point. "
                "The entry (i,j) is partial_s^i partial_t^j K(s0,s0)/(i!j!). "
                "If the kernel were positive, every such matrix would be positive semidefinite."
            ),
            "lambda": (
                "The lowest eigenvector of the 9-jet matrix gives a local linear test "
                "Lambda_a(f)=sum e_k(a) f^(k)(a)/k!.  The active endpoint obstruction "
                "is the moving family of these tests."
            ),
            "schur": (
                "R collects the Lambda_a traces.  Splitting functions into ker R and its "
                "orthogonal complement turns the proof into a Schur/Douglas range problem."
            ),
        },
    }

    Path("visual_data.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    print("wrote visual_data.json")


if __name__ == "__main__":
    main()
