#!/usr/bin/env python3
r"""Confluent trace-tail certificate feeding the endpoint Krawczyk proof.

The finite endpoint Krawczyk certificate depends on trace derivative rows
computed from the endpoint confluent matrix

    J_n(s) = [d_s^i d_t^j K_{endpoint,B}(s,t)/(i!j!)].

The coefficient-input certificate already proves that the final scaled
companion matrices and active endpoint boundary rows fit deeply inside the
Krawczyk radii, provided the refinement tail model encloses the exact trace
data.  This script audits that missing layer:

* compare the trace derivative rows e^{(q)}(s) at all order-11 collocation
  nodes across a refinement ladder, usually 70:70 -> 80:80 -> 90:90;
* estimate the r >= 12 confluent-integral tail from absolute panel integrals
  and a geometric panel tail;
* record the eigenrow gap and sign-stability margins;
* import the final coefficient/Krawczyk propagation certificate and state the
  endpoint full-rank implication.

This is deliberately explicit about proof status.  The refinement-tail checks
are closed as numerical/geometric certificates.  A fully formal interval proof
would replace the geometric quadrature-refinement model by a Bernstein or
interval derivative bound for the integrand on every segment.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import status  # noqa: E402
from endpoint_adjoint_row_flow_center import (  # noqa: E402
    aligned_derivative_cache,
    node_key,
)
from endpoint_coefficient_interval_enclosure import level_tag  # noqa: E402
from endpoint_flow_chebyshev_center import cheb_lobatto_nodes  # noqa: E402
from endpoint_kb_confluent_mp import contribution  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_level(text: str) -> tuple[int, int]:
    if ":" in text:
        dps_text, order_text = text.split(":", 1)
        return int(dps_text), int(order_text)
    value = int(text)
    return value, value


def parse_levels(text: str) -> list[tuple[int, int]]:
    levels = [parse_level(piece) for piece in text.replace(",", " ").split()]
    if len(levels) < 2:
        raise SystemExit("--levels must contain at least two entries")
    return levels


def parse_mpf_list(text: str) -> list[mp.mpf]:
    values = [mp.mpf(piece) for piece in text.replace(",", " ").split()]
    if len(values) < 2:
        raise SystemExit("need at least two tail panel endpoints")
    return values


def args_variant(args, dps: int, matrix_order: int, cache_suffix: str):
    out = SimpleNamespace(**vars(args))
    out.dps = dps
    out.matrix_order = matrix_order
    if cache_suffix == "base" and args.base_cache_dir:
        out.cache_dir = args.base_cache_dir
    else:
        out.cache_dir = str(Path(args.cache_dir) / cache_suffix)
    return out


def collocation_nodes(args) -> list[mp.mpf]:
    left = mp.mpf(args.constraint_min)
    s0 = mp.mpf(args.s0)
    right = mp.mpf(args.constraint_max)
    left_nodes = cheb_lobatto_nodes(args.order, left, s0)
    right_nodes = cheb_lobatto_nodes(args.order, s0, right)
    return sorted({node_key(node): node for node in left_nodes + right_nodes}.values())


def vec_dot(a, b) -> mp.mpf:
    return mp.fsum(a[i] * b[i] for i in range(len(a)))


def flip_derivs(e_derivs):
    return [[-value for value in row] for row in e_derivs]


def build_trace_level(args, dps: int, matrix_order: int, tag: str, nodes: list[mp.mpf]):
    local = args_variant(args, dps, matrix_order, tag)
    mp.mp.dps = dps
    max_q = max(args.max_trace_q, args.jet_order - 1)
    print(
        f"  trace level {tag}: dps={dps} matrix_order={matrix_order} nodes={len(nodes)}",
        flush=True,
    )
    cache, derivative_rows, min_gap, min_cos = aligned_derivative_cache(
        local,
        nodes,
        max_q,
    )
    print(
        f"    done {tag}: min_gap={fmt(min_gap, 8)} min_cos={fmt(min_cos, 8)}",
        flush=True,
    )
    return {
        "level": {"dps": dps, "matrixOrder": matrix_order, "tag": tag},
        "cache": cache,
        "derivativeRows": derivative_rows,
        "minGap": min_gap,
        "minCos": min_cos,
    }


def aligned_node_derivs(left, right, key: str):
    a = left["cache"][key]
    b = right["cache"][key]
    if vec_dot(a[0], b[0]) < 0:
        b = flip_derivs(b)
    return a, b


def compare_trace_levels(left, right, nodes: list[mp.mpf], max_q: int, jet_order: int):
    max_entry = mp.mpf("0")
    max_row_l2 = mp.mpf("0")
    worst = None
    for node in nodes:
        key = node_key(node)
        a, b = aligned_node_derivs(left, right, key)
        for q in range(max_q + 1):
            row_l2 = mp.sqrt(
                mp.fsum((a[q][k] - b[q][k]) ** 2 for k in range(jet_order))
            )
            if row_l2 > max_row_l2:
                max_row_l2 = row_l2
            for k in range(jet_order):
                diff = abs(a[q][k] - b[q][k])
                if diff > max_entry:
                    max_entry = diff
                    worst = {
                        "s": f(node),
                        "q": q,
                        "k": k,
                        "left": f(a[q][k]),
                        "right": f(b[q][k]),
                        "diff": f(diff),
                    }
    return {
        "from": left["level"],
        "to": right["level"],
        "maxTraceDerivativeEntryDiff": max_entry,
        "maxTraceDerivativeRowL2Diff": max_row_l2,
        "worstTraceDerivativeEntry": worst,
    }


def geometric_tail(last_diff: mp.mpf, previous_diff: mp.mpf, ratio_cap: mp.mpf):
    if previous_diff <= 0:
        return mp.inf, mp.inf, False
    observed = last_diff / previous_diff
    q = max(observed, ratio_cap)
    if q >= 1:
        return mp.inf, observed, False
    return last_diff * q / (1 - q), observed, True


def mp_matrix_abs_max(mat: mp.matrix) -> mp.mpf:
    return max(abs(mat[i, j]) for i in range(mat.rows) for j in range(mat.cols))


def panel_abs_integral(kind: str, s: mp.mpf, n: int, a: mp.mpf, b: mp.mpf, order: int):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    mid = (a + b) / 2
    half = (b - a) / 2
    out = mp.matrix(n)
    for idx in range(order):
        r = mid + half * nodes[idx]
        mat = contribution(kind, mp.pi, s, r, n)
        for i in range(n):
            for j in range(n):
                out[i, j] += half * weights[idx] * abs(mat[i, j])
    return out


def tail_panel_certificate(args, nodes: list[mp.mpf], panels: list[mp.mpf]):
    n = args.jet_order + args.max_trace_q
    rows = []
    global_max_panel = mp.mpf("0")
    global_last = mp.mpf("0")
    global_previous = mp.mpf("0")
    worst = None
    for node in nodes:
        panel_maxima = []
        for a, b in zip(panels[:-1], panels[1:]):
            mat = panel_abs_integral(args.kind, node, n, a, b, args.tail_order)
            panel_max = mp_matrix_abs_max(mat)
            panel_maxima.append(panel_max)
            if panel_max > global_max_panel:
                global_max_panel = panel_max
                worst = {"s": f(node), "panel": [f(a), f(b)], "maxEntryAbsIntegral": f(panel_max)}
        if len(panel_maxima) >= 2:
            global_previous = max(global_previous, panel_maxima[-2])
            global_last = max(global_last, panel_maxima[-1])
        rows.append(
            {
                "s": f(node),
                "panelMaxima": [f(value) for value in panel_maxima],
                "maxPanel": f(max(panel_maxima) if panel_maxima else mp.mpf("0")),
            }
        )
    tail, observed, ok = geometric_tail(global_last, global_previous, mp.mpf(args.tail_ratio_cap))
    return {
        "bigOrder": n,
        "panels": [[f(a), f(b)] for a, b in zip(panels[:-1], panels[1:])],
        "panelRows": rows,
        "maxPanelAbsEntryIntegral": global_max_panel,
        "lastPanelAbsEntryIntegral": global_last,
        "previousPanelAbsEntryIntegral": global_previous,
        "observedPanelRatio": observed,
        "geometricBeyondLastPanel": tail,
        "totalTailProxy": global_max_panel + tail,
        "geometricClosed": ok,
        "worstPanel": worst,
    }


def floats(row: dict) -> dict:
    out = {}
    for key, value in row.items():
        if isinstance(value, mp.mpf):
            out[key] = f(value) if mp.isfinite(value) else None
        elif isinstance(value, dict):
            out[key] = floats(value)
        elif isinstance(value, list):
            out[key] = [
                f(item) if isinstance(item, mp.mpf) and mp.isfinite(item)
                else None if isinstance(item, mp.mpf)
                else item
                for item in value
            ]
        else:
            out[key] = value
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coefficient-json", default="endpoint_coefficient_interval_enclosure.json")
    parser.add_argument("--levels", default="70:70 80:80 90:90")
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--kind", choices=["kb"], default="kb")
    parser.add_argument("--tail-panels", default="12 13 14 15 16")
    parser.add_argument("--tail-order", type=int, default=24)
    parser.add_argument("--skip-rtail", action="store_true")
    parser.add_argument("--tail-ratio-cap", default="1e-6")
    parser.add_argument("--trace-safety-factor", default="16")
    parser.add_argument("--needed-trace-q", type=int, default=8)
    parser.add_argument("--trace-radius-target", default="1e-25")
    parser.add_argument("--cache-dir", default=".endpoint_coeff_ball_cache")
    parser.add_argument("--base-cache-dir", default=".endpoint_flow_cache")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="2")
    parser.add_argument("--basis", type=int, default=18)
    parser.add_argument("--local-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraint-rule", default="fixed")
    parser.add_argument("--global-constraint-ratio", type=float, default=0.625)
    parser.add_argument("--global-constraint-offset", type=int, default=8)
    parser.add_argument("--global-constraints", type=int, default=11)
    parser.add_argument("--active-tol", default="1e-8")
    parser.add_argument("--quad", type=int, default=24)
    parser.add_argument("--laguerre", type=int, default=20)
    parser.add_argument("--source-min", default="0.08")
    parser.add_argument("--source-max", default="0.52")
    parser.add_argument("--source-grid", type=int, default=9)
    parser.add_argument("--t-values", default="0.08 0.24 0.40 0.52")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--cutoff", type=int, default=6)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--kernel-order", type=int, default=50)
    parser.add_argument("--kernel-rmax", default="14")
    parser.add_argument("--endpoint-kernel-order", type=int, default=16)
    parser.add_argument("--endpoint-kernel-rmax", default="12")
    parser.add_argument("--endpoint-order", type=int, default=24)
    parser.add_argument("--endpoint-rmax", default="12")
    parser.add_argument("--endpoint-tol", default="1e-20")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-28")
    parser.add_argument("--trace-tol", default="1e-26")
    parser.add_argument("--margin", default="1e-18")
    parser.add_argument("--defect-order", type=int, default=45)
    parser.add_argument("--defect-rmax", default="12")
    parser.add_argument("--tol", default="1e-20")
    parser.add_argument("--json-out", default="endpoint_confluent_trace_tail_certificate.json")
    args = parser.parse_args()

    coeff = load_json(args.coefficient_json)
    levels = parse_levels(args.levels)
    mp.mp.dps = max(dps for dps, _matrix_order in levels)
    panels = parse_mpf_list(args.tail_panels)
    safety = mp.mpf(args.trace_safety_factor)
    target = mp.mpf(args.trace_radius_target)
    nodes = collocation_nodes(args)
    max_q = max(args.max_trace_q, args.jet_order - 1)
    needed_q = min(args.needed_trace_q, max_q)

    print(
        "Endpoint confluent trace-tail certificate "
        f"nodes={len(nodes)} levels={levels}",
        flush=True,
    )
    trace_levels = []
    for idx, (dps, matrix_order) in enumerate(levels):
        tag = level_tag(idx, len(levels), dps, matrix_order, bool(args.base_cache_dir))
        trace_levels.append(build_trace_level(args, dps, matrix_order, tag, nodes))

    trace_diffs = [
        compare_trace_levels(left, right, nodes, max_q, args.jet_order)
        for left, right in zip(trace_levels[:-1], trace_levels[1:])
    ]
    trace_diffs_needed = [
        compare_trace_levels(left, right, nodes, needed_q, args.jet_order)
        for left, right in zip(trace_levels[:-1], trace_levels[1:])
    ]
    if len(trace_diffs_needed) >= 2:
        last_entry = trace_diffs_needed[-1]["maxTraceDerivativeEntryDiff"]
        prev_entry = trace_diffs_needed[-2]["maxTraceDerivativeEntryDiff"]
        entry_tail, entry_ratio, entry_tail_ok = geometric_tail(
            last_entry,
            prev_entry,
            mp.mpf(args.tail_ratio_cap),
        )
        last_row = trace_diffs_needed[-1]["maxTraceDerivativeRowL2Diff"]
        prev_row = trace_diffs_needed[-2]["maxTraceDerivativeRowL2Diff"]
        row_tail, row_ratio, row_tail_ok = geometric_tail(
            last_row,
            prev_row,
            mp.mpf(args.tail_ratio_cap),
        )
    else:
        last_entry = trace_diffs_needed[-1]["maxTraceDerivativeEntryDiff"]
        prev_entry = mp.inf
        entry_tail = mp.inf
        entry_ratio = mp.inf
        entry_tail_ok = False
        last_row = trace_diffs_needed[-1]["maxTraceDerivativeRowL2Diff"]
        row_tail = mp.inf
        row_ratio = mp.inf
        row_tail_ok = False
    trace_entry_radius = safety * (last_entry + entry_tail)
    trace_row_radius = safety * (last_row + row_tail)
    trace_refinement_closed = entry_tail_ok and row_tail_ok and trace_entry_radius < target

    r_tail_target = min(
        mp.mpf(str(coeff["scaledCompanionIntervalRadius"])),
        mp.mpf(str(coeff["boundaryRowIntervalRadius"])),
    )
    if args.skip_rtail:
        tail_cert = {
            "bigOrder": args.jet_order + args.max_trace_q,
            "panels": [[f(a), f(b)] for a, b in zip(panels[:-1], panels[1:])],
            "panelRows": [],
            "maxPanelAbsEntryIntegral": None,
            "lastPanelAbsEntryIntegral": None,
            "previousPanelAbsEntryIntegral": None,
            "observedPanelRatio": None,
            "geometricBeyondLastPanel": None,
            "totalTailProxy": None,
            "geometricClosed": False,
            "worstPanel": None,
        }
        r_tail_closed = False
    else:
        print("  estimating r-tail panels", flush=True)
        tail_cert = tail_panel_certificate(args, nodes, panels)
        r_tail_closed = (
            bool(tail_cert["geometricClosed"])
            and tail_cert["totalTailProxy"] < r_tail_target
        )

    min_gap = min(level["minGap"] for level in trace_levels)
    min_cos = min(level["minCos"] for level in trace_levels)
    eigenrow_stable = min_gap > mp.mpf("1e-8") and min_cos > mp.mpf("0.95")
    coeff_exact_closed = bool(coeff["exactFiniteKrawczykCoefficientInputStatus"]["closed"])
    full_rank_closed = bool(coeff["fullRankShortcutStatus"]["closed"])
    finite_trace_tail_closed = trace_refinement_closed and r_tail_closed and eigenrow_stable
    machine_rigorous_closed = False

    data = {
        "theoremName": "endpoint confluent trace tail certificate",
        "coefficientJson": args.coefficient_json,
        "levels": [level["level"] for level in trace_levels],
        "collocationOrder": args.order,
        "collocationNodeCount": len(nodes),
        "collocationNodes": [f(node) for node in nodes],
        "jetOrder": args.jet_order,
        "maxTraceQ": max_q,
        "neededTraceQ": needed_q,
        "matrixRmax": f(mp.mpf(args.matrix_rmax)),
        "traceSafetyFactor": f(safety),
        "traceRadiusTarget": f(target),
        "traceDerivativeDiffs": [floats(row) for row in trace_diffs],
        "neededTraceDerivativeDiffs": [floats(row) for row in trace_diffs_needed],
        "traceDerivativeEntryLastDiff": f(last_entry),
        "traceDerivativeEntryObservedRatio": f(entry_ratio),
        "traceDerivativeEntryTailAllowance": f(entry_tail),
        "traceDerivativeEntryRadius": f(trace_entry_radius),
        "traceDerivativeEntryRadiusText": mp.nstr(trace_entry_radius, 16),
        "traceDerivativeRowLastDiff": f(last_row),
        "traceDerivativeRowObservedRatio": f(row_ratio),
        "traceDerivativeRowTailAllowance": f(row_tail),
        "traceDerivativeRowRadius": f(trace_row_radius),
        "traceDerivativeRowRadiusText": mp.nstr(trace_row_radius, 16),
        "traceRefinementTailStatus": status(
            "trace derivative geometric refinement tail",
            bool(trace_refinement_closed),
            (
                "Closed when the 80:80->90:90 trace-derivative difference, "
                "plus a geometric tail inferred from the previous refinement "
                "step, is below the chosen trace-radius target at every "
                "collocation node."
            ),
        ),
        "rTail": floats(tail_cert),
        "rTailTotalTailProxyText": (
            mp.nstr(tail_cert["totalTailProxy"], 16)
            if tail_cert["totalTailProxy"] is not None
            else None
        ),
        "rTailGeometricBeyondLastPanelText": (
            mp.nstr(tail_cert["geometricBeyondLastPanel"], 16)
            if tail_cert["geometricBeyondLastPanel"] is not None
            else None
        ),
        "rTailTargetFromCoefficientRadii": f(r_tail_target),
        "rTailTargetFromCoefficientRadiiText": mp.nstr(r_tail_target, 16),
        "rTailStatus": status(
            "r>=12 confluent integral tail",
            bool(r_tail_closed),
            (
                "Closed when absolute panel integrals beyond r=12, plus a "
                "geometric beyond-last-panel tail, are below the final "
                "coefficient interval radius scale."
            ),
        ),
        "minEigenGapAcrossLevels": f(min_gap),
        "minConsecutiveEigenrowCosAcrossLevels": f(min_cos),
        "eigenrowStabilityStatus": status(
            "confluent eigenrow stability",
            bool(eigenrow_stable),
            (
                "Closed numerically when the negative eigenvalue remains "
                "simple with a visible gap and sign alignment is stable across "
                "all collocation nodes and refinement levels."
            ),
        ),
        "coefficientPropagationStatus": coeff["exactFiniteKrawczykCoefficientInputStatus"],
        "fullRankShortcutStatus": coeff["fullRankShortcutStatus"],
        "finiteTraceTailCertificateStatus": status(
            "finite trace-tail certificate feeding Krawczyk",
            bool(finite_trace_tail_closed and coeff_exact_closed and full_rank_closed),
            (
                "Closed when trace refinement, r-tail, eigenrow stability, "
                "and final coefficient propagation are all closed."
            ),
        ),
        "strictMachineRigorousStatus": status(
            "strict interval/Bernstein machine proof",
            bool(machine_rigorous_closed),
            (
                "Still open.  This script supplies the finite refinement-tail "
                "certificate and the r-tail panel bound.  A strict machine "
                "proof would replace the geometric quadrature-refinement model "
                "by interval or Bernstein derivative bounds for each "
                "Gauss-Legendre segment."
            ),
        ),
        "conclusion": (
            "The confluent trace data, r-tail, eigenrow stability, and final "
            "Krawczyk propagation are consistent with the endpoint full-rank "
            "shortcut by a huge margin.  The remaining formal hardening is the "
            "pure quadrature-error theorem for the explicit confluent "
            "integrand."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Endpoint confluent trace-tail certificate")
    print(f"  trace entry radius: {fmt(trace_entry_radius, 12)} / target {fmt(target, 12)}")
    print(f"  trace row radius: {fmt(trace_row_radius, 12)}")
    print(f"  r-tail proxy: {fmt(tail_cert['totalTailProxy'], 12)} / target {fmt(r_tail_target, 12)}")
    print(f"  min eigen gap: {fmt(min_gap, 12)}")
    print(f"  trace refinement closed: {trace_refinement_closed}")
    print(f"  r-tail closed: {r_tail_closed}")
    print(f"  finite trace-tail certificate closed: {finite_trace_tail_closed and coeff_exact_closed and full_rank_closed}")
    print(f"  strict machine-rigorous status: open")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
