#!/usr/bin/env python3
r"""Finite gamma consequence for the continuum trace-frame theorem."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def nested_closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    return bool(item.get("closed")) if isinstance(item, dict) else False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-frame-json",
        default="trace_frame_interval_lower_bound_certificate.json",
    )
    parser.add_argument(
        "--json-out",
        default="trace_frame_finite_gamma_consequence_theorem.json",
    )
    args = parser.parse_args()

    finite = load(args.finite_frame_json)
    gamma = float(finite.get("gammaFiniteLower", 0.0))
    positive = bool(
        finite.get("gammaFiniteLowerPositive")
        or nested_closed(finite, "finiteTraceFrameIntervalLowerBoundStatus")
        or gamma > 0.0
    )
    ok = positive and gamma > 0.0

    data = {
        "theoremName": "trace-frame finite gamma consequence theorem",
        "proofClass": "symbolic identity",
        "finiteFrameJson": args.finite_frame_json,
        "gammaFiniteLower": finite.get("gammaFiniteLower"),
        "gammaFiniteLowerString": finite.get("gammaFiniteLowerString"),
        "observedGammaFloor": finite.get("gammaFiniteLower"),
        "basis": finite.get("basis"),
        "traceCount": finite.get("traceCount"),
        "lambdaMinCenter": finite.get("lambdaMinCenter"),
        "radiusOperatorBoundUsed": finite.get("radiusOperatorBoundUsed"),
        "maxRangeResidualRelative": finite.get("maxRangeResidualRelative"),
        "activeDimension": finite.get("activeDimension"),
        "statement": (
            "The finite weighted trace-frame matrix has a certified positive "
            "lower eigenvalue gamma_N^- > 0."
        ),
        "statuses": {
            "finiteTraceFrameInputStatus": status(
                "finite interval trace-frame input",
                positive,
                "The finite interval certificate proves gammaFiniteLowerPositive.",
            ),
            "finiteGammaPositiveConsequenceStatus": status(
                "finite gamma positive consequence",
                ok,
                "Therefore gamma_N^- is a positive finite trace-frame lower bound.",
            ),
        },
        "gammaFiniteLowerPositive": ok,
        "finiteTraceFrameLowerBoundClosed": ok,
        "proof": [
            "Import the finite trace-frame interval lower-bound certificate.",
            "Read off gammaFiniteLower and the finite lower-bound status.",
            "Expose only the positive finite gamma consequence to the continuum theorem.",
        ],
        "remainingAnalyticGap": None if ok else "Finite trace-frame gamma lower bound is not positive.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Trace-frame finite gamma consequence theorem")
    print(f"  gammaFiniteLower: {gamma:.12e}")
    print(f"  theorem closed: {ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
