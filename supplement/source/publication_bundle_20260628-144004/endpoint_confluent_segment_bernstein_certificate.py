#!/usr/bin/env python3
r"""Bernstein segment bounds for the endpoint confluent Gauss quadrature.

This is the hardening layer beneath
``endpoint_confluent_trace_tail_certificate.py``.  The previous certificate
used a refinement ladder to model the finite Gauss-Legendre error.  Here we
bound the error on explicit Gauss-Legendre segments.

For a segment ``r = mid + half*t`` and an analytic matrix entry ``g(t)``,
write its Chebyshev expansion on ``[-1,1]``.  If ``|g| <= M_rho`` on the
Bernstein ellipse ``E_rho``, then the Chebyshev coefficients satisfy
``|c_k| <= 2 M_rho rho^{-k}``.  Since N-point Gauss-Legendre quadrature is
exact for all polynomials of degree <= 2N-1,

    | int g - Q_N g |
      <= half * 8 M_rho rho^{-2N} / (1-rho^{-1}).

The bound is intentionally conservative: ``|int T_k - Q_N T_k| <= 4`` for
the tail.  The payoff is that it is simple and local to each explicit segment.

The default ``majorant`` evaluator does not use dependency-heavy interval
arithmetic on the whole formula.  It propagates absolute majorants through the
exact confluent Taylor recurrence, retaining the lower bound on
``Re(exp(r)-1)`` that supplies the large-r exponential decay.  This is the
deterministic Bernstein replacement for the geometric quadrature-refinement
model.  The older ``interval`` and ``ball`` evaluators remain available as
diagnostics.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from adjoint_green_endpoint_selection import status  # noqa: E402
from endpoint_confluent_trace_tail_certificate import collocation_nodes, parse_mpf_list  # noqa: E402
from endpoint_kb_confluent_mp import segments  # noqa: E402
from global_trace_observability_gap import f, fmt  # noqa: E402


class CBall:
    __slots__ = ("center", "radius")

    def __init__(self, center=0, radius=0):
        self.center = mp.mpc(center)
        self.radius = mp.mpf(radius)

    @staticmethod
    def as_ball(value):
        if isinstance(value, CBall):
            return value
        return CBall(value, 0)

    def __add__(self, other):
        other = CBall.as_ball(other)
        return CBall(self.center + other.center, self.radius + other.radius)

    __radd__ = __add__

    def __sub__(self, other):
        other = CBall.as_ball(other)
        return CBall(self.center - other.center, self.radius + other.radius)

    def __rsub__(self, other):
        other = CBall.as_ball(other)
        return CBall(other.center - self.center, self.radius + other.radius)

    def __neg__(self):
        return CBall(-self.center, self.radius)

    def __mul__(self, other):
        other = CBall.as_ball(other)
        radius = (
            abs(self.center) * other.radius
            + abs(other.center) * self.radius
            + self.radius * other.radius
        )
        return CBall(self.center * other.center, radius)

    __rmul__ = __mul__

    def inv(self):
        mag = abs(self.center)
        if mag <= self.radius:
            return CBall(mp.inf, mp.inf)
        return CBall(1 / self.center, self.radius / (mag * (mag - self.radius)))

    def __truediv__(self, other):
        other = CBall.as_ball(other)
        return self * other.inv()

    def __rtruediv__(self, other):
        other = CBall.as_ball(other)
        return other * self.inv()

    def exp(self):
        exp_center = mp.e ** self.center
        radius = mp.e ** (mp.re(self.center)) * (mp.e ** self.radius - 1)
        return CBall(exp_center, radius)

    def abs_sup(self):
        return abs(self.center) + self.radius


def load_json(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def iv_upper(x) -> mp.mpf:
    return mp.mpf(x.b)


def iv_lower(x) -> mp.mpf:
    return mp.mpf(x.a)


def interval_abs_sup(x) -> mp.mpf:
    if isinstance(x, CBall):
        return x.abs_sup()
    name = type(x).__name__
    if hasattr(x, "_mpci_") or name.endswith("ivmpc"):
        return iv_upper(abs(x))
    if hasattr(x, "_mpi_") or name.endswith("ivmpf"):
        lo = iv_lower(x)
        hi = iv_upper(x)
        return max(abs(lo), abs(hi))
    return abs(x)


def ellipse_arc_box(theta_left, theta_right, rho):
    theta = mp.iv.mpf([theta_left, theta_right])
    rho_iv = mp.iv.mpf([rho, rho])
    major = (rho_iv + 1 / rho_iv) / 2
    minor = (rho_iv - 1 / rho_iv) / 2
    return mp.iv.mpc(major * mp.iv.cos(theta), minor * mp.iv.sin(theta))


def ellipse_arc_rect_box(theta_left, theta_right, rho):
    """Rectangular enclosure of a Bernstein ellipse arc.

    Direct interval evaluation of ``cos([a,b])`` and ``sin([a,b])`` creates
    severe dependency blowup inside the endpoint Taylor series.  This ball
    enclosure uses the arc midpoint and the elementary Lipschitz bound
    ``|z'(theta)| <= (rho+rho^{-1})/2``.
    """
    mid_theta = (theta_left + theta_right) / 2
    half_width = (theta_right - theta_left) / 2
    major = (rho + 1 / rho) / 2
    minor = (rho - 1 / rho) / 2
    center_real = major * mp.cos(mid_theta)
    center_imag = minor * mp.sin(mid_theta)
    radius = major * half_width
    real = mp.iv.mpf([center_real - radius, center_real + radius])
    imag = mp.iv.mpf([center_imag - radius, center_imag + radius])
    return mp.iv.mpc(real, imag)


def ellipse_arc_disk(theta_left, theta_right, rho):
    mid_theta = (theta_left + theta_right) / 2
    half_width = (theta_right - theta_left) / 2
    major = (rho + 1 / rho) / 2
    minor = (rho - 1 / rho) / 2
    center = major * mp.cos(mid_theta) + 1j * minor * mp.sin(mid_theta)
    radius = major * half_width
    return center, radius


def rho_from_imag_limit(half: mp.mpf, max_imag: mp.mpf, rho_cap: mp.mpf) -> mp.mpf:
    if half <= 0:
        return rho_cap
    d = 2 * max_imag / half
    rho = (d + mp.sqrt(d * d + 4)) / 2
    return min(rho, rho_cap)


def segment_error_bound(length: mp.mpf, rho: mp.mpf, order: int, mrho: mp.mpf) -> mp.mpf:
    half = length / 2
    return half * 8 * mrho * rho ** (-2 * order) / (1 - 1 / rho)


def refined_segments(rmax: mp.mpf, step: mp.mpf) -> list[mp.mpf]:
    if step <= 0:
        return segments(rmax)
    pts = [mp.mpf("0")]
    cur = mp.mpf("0")
    guard = 0
    while cur < rmax:
        nxt = min(rmax, cur + step)
        if nxt <= cur:
            raise ValueError("segment step did not advance")
        pts.append(nxt)
        cur = nxt
        guard += 1
        if guard > 100000:
            raise ValueError("too many refined segments")
    return pts


def sum_any(values):
    total = None
    for value in values:
        total = value if total is None else total + value
    return mp.mpf("0") if total is None else total


def zero(n: int):
    return [mp.mpf("0") for _ in range(n)]


def zero_ball(n: int):
    return [CBall(0, 0) for _ in range(n)]


def sub_poly(a, b):
    return [x - y for x, y in zip(a, b)]


def scale_poly(a, c):
    return [c * x for x in a]


def mul_poly(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = sum_any(a[j] * b[i - j] for j in range(i + 1))
    return out


def div_poly(a, b):
    n = len(a)
    out = zero(n)
    for i in range(n):
        out[i] = (
            a[i] - sum_any(b[j] * out[i - j] for j in range(1, i + 1))
        ) / b[0]
    return out


def exp_series_any(a):
    n = len(a)
    out = zero(n)
    out[0] = mp.e ** a[0]
    for i in range(1, n):
        out[i] = sum_any(k * a[k] * out[i - k] for k in range(1, i + 1)) / i
    return out


def exp_series_ball(a):
    n = len(a)
    out = zero_ball(n)
    out[0] = CBall.as_ball(a[0]).exp()
    for i in range(1, n):
        out[i] = sum_any(k * a[k] * out[i - k] for k in range(1, i + 1)) / i
    return out


def lambda_series(c, s0, n):
    lam0 = c * mp.e ** s0
    return [lam0 / mp.factorial(i) for i in range(n)]


def endpoint_series_any(c, s0, r, n):
    x = mp.e ** r
    tau = x - 1
    lam = lambda_series(c, s0, n)
    denom = lam[:]
    denom[0] -= c

    exp_lam = exp_series_any(scale_poly(lam, -tau))
    exp_c = mp.e ** (-c * tau)

    one = zero(n)
    one[0] = 1
    sser = zero(n)
    sser[0] = s0
    if n > 1:
        sser[1] = 1

    f_num = exp_lam[:]
    f_num[0] -= exp_c
    f = div_poly(f_num, denom)

    gb_pref = sub_poly(one, scale_poly(lam, x))
    gb = mul_poly(gb_pref, exp_lam)
    gb_c = (1 - c * x) * exp_c
    vb_num = gb[:]
    vb_num[0] -= gb_c
    vb = div_poly(vb_num, denom)
    wb = div_poly(mul_poly(sser, gb), denom)

    return f, vb, wb


def endpoint_series_ball(c, s0, r, n):
    x = r.exp()
    tau = x - 1
    lam = [CBall(value, 0) for value in lambda_series(c, s0, n)]
    denom = lam[:]
    denom[0] = denom[0] - c

    exp_lam = exp_series_ball(scale_poly(lam, -tau))
    exp_c = (-c * tau).exp()

    one = zero_ball(n)
    one[0] = CBall(1, 0)
    sser = zero_ball(n)
    sser[0] = CBall(s0, 0)
    if n > 1:
        sser[1] = CBall(1, 0)

    f_num = exp_lam[:]
    f_num[0] = f_num[0] - exp_c
    f = div_poly(f_num, denom)

    gb_pref = sub_poly(one, scale_poly(lam, x))
    gb = mul_poly(gb_pref, exp_lam)
    gb_c = (1 - c * x) * exp_c
    vb_num = gb[:]
    vb_num[0] = vb_num[0] - gb_c
    vb = div_poly(vb_num, denom)
    wb = div_poly(mul_poly(sser, gb), denom)
    return f, vb, wb


def contribution_any(kind, c, s0, r, n):
    if kind != "kb":
        raise ValueError(kind)
    _f, vb, wb = endpoint_series_any(c, s0, r, n)
    mat = [[mp.mpf("0") for _ in range(n)] for _ in range(n)]
    weight = mp.e ** (mp.mpf("2.5") * r)
    for i in range(n):
        for j in range(i + 1):
            val = r * vb[i] * vb[j] + mp.mpf("0.5") * (
                vb[i] * wb[j] + wb[i] * vb[j]
            )
            mat[i][j] = mat[j][i] = weight * val
    return mat


def contribution_ball(kind, c, s0, r, n):
    if kind != "kb":
        raise ValueError(kind)
    _f, vb, wb = endpoint_series_ball(c, s0, r, n)
    mat = [[CBall(0, 0) for _ in range(n)] for _ in range(n)]
    weight = (mp.mpf("2.5") * r).exp()
    for i in range(n):
        for j in range(i + 1):
            val = r * vb[i] * vb[j] + mp.mpf("0.5") * (
                vb[i] * wb[j] + wb[i] * vb[j]
            )
            mat[i][j] = mat[j][i] = weight * val
    return mat


def div_bound(num, denom_tail, denom0):
    n = len(num)
    out = [mp.mpf("0") for _ in range(n)]
    for i in range(n):
        out[i] = (
            num[i] + mp.fsum(denom_tail[j] * out[i - j] for j in range(1, i + 1))
        ) / denom0
    return out


def endpoint_series_majorants(c, s0, a, b, rho, n):
    """Absolute majorants for the B-branch endpoint Taylor coefficients.

    The Taylor variable is the trace/base variable used by
    ``endpoint_kb_confluent_mp.endpoint_series``.  The complex variable being
    enclosed is the segment coordinate on the Bernstein ellipse.
    """
    mid = (a + b) / 2
    half = (b - a) / 2
    major = (rho + 1 / rho) / 2
    minor = (rho - 1 / rho) / 2
    re_r_max = mid + half * major
    re_r_min = mid - half * major
    im_r_max = abs(half) * minor
    x_abs = mp.e ** re_r_max
    tau_abs = x_abs + 1
    re_tau_min = mp.e ** re_r_min * mp.cos(im_r_max) - 1

    lam = [c * mp.e ** s0 / mp.factorial(k) for k in range(n)]
    denom0 = lam[0] - c
    if denom0 <= 0:
        raise ValueError("endpoint denominator is not bounded away from zero")
    denom_tail = [abs(value) for value in lam]

    exp_lam = [mp.mpf("0") for _ in range(n)]
    exp_lam[0] = mp.e ** (-lam[0] * re_tau_min)
    for i in range(1, n):
        exp_lam[i] = mp.fsum(
            k * tau_abs * abs(lam[k]) * exp_lam[i - k]
            for k in range(1, i + 1)
        ) / i
    exp_c = mp.e ** (-c * re_tau_min)

    gb_pref = [mp.mpf("0") for _ in range(n)]
    gb_pref[0] = 1 + abs(lam[0]) * x_abs
    for i in range(1, n):
        gb_pref[i] = abs(lam[i]) * x_abs
    gb = [mp.mpf("0") for _ in range(n)]
    for i in range(n):
        gb[i] = mp.fsum(gb_pref[j] * exp_lam[i - j] for j in range(i + 1))

    vb_num = gb[:]
    vb_num[0] += (1 + c * x_abs) * exp_c
    vb = div_bound(vb_num, denom_tail, denom0)

    sser = [mp.mpf("0") for _ in range(n)]
    sser[0] = abs(s0)
    if n > 1:
        sser[1] = 1
    wb_num = [mp.mpf("0") for _ in range(n)]
    for i in range(n):
        wb_num[i] = mp.fsum(sser[j] * gb[i - j] for j in range(min(i + 1, 2)))
    wb = div_bound(wb_num, denom_tail, denom0)

    r_abs = abs(mid) + abs(half) * major
    weight_abs = mp.e ** (mp.mpf("2.5") * re_r_max)
    return vb, wb, r_abs, weight_abs, re_tau_min


def contribution_majorant(kind, c, s0, a, b, rho, n):
    if kind != "kb":
        raise ValueError(kind)
    vb, wb, r_abs, weight_abs, re_tau_min = endpoint_series_majorants(c, s0, a, b, rho, n)
    mat = [[mp.mpf("0") for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            val = r_abs * vb[i] * vb[j] + mp.mpf("0.5") * (
                vb[i] * wb[j] + wb[i] * vb[j]
            )
            mat[i][j] = mat[j][i] = weight_abs * val
    return mat, re_tau_min


def matrix_abs_bound(mat) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    entry = mp.mpf("0")
    row = mp.mpf("0")
    frob2 = mp.mpf("0")
    row_count = mat.rows if hasattr(mat, "rows") else len(mat)
    col_count = mat.cols if hasattr(mat, "cols") else len(mat[0])
    for i in range(row_count):
        row_sum = mp.mpf("0")
        for j in range(col_count):
            value = mat[i, j] if hasattr(mat, "rows") else mat[i][j]
            bound = interval_abs_sup(value)
            entry = max(entry, bound)
            row_sum += bound
            frob2 += bound * bound
        row = max(row, row_sum)
    return entry, row, mp.sqrt(frob2)


def segment_node_mrho(args, s: mp.mpf, a: mp.mpf, b: mp.mpf, rho: mp.mpf):
    mid = (a + b) / 2
    half = (b - a) / 2
    n = args.jet_order + args.needed_trace_q
    if args.eval_mode == "majorant":
        mat, re_tau_min = contribution_majorant(args.kind, mp.pi, s, a, b, rho, n)
        entry, row, frob = matrix_abs_bound(mat)
        return entry, row, frob, {
            "type": "absolute-majorant",
            "reTauMin": f(re_tau_min),
        }
    max_entry = mp.mpf("0")
    max_row = mp.mpf("0")
    max_frob = mp.mpf("0")
    worst_arc = None
    step = 2 * mp.pi / args.arcs
    for arc in range(args.arcs):
        theta_left = arc * step
        theta_right = (arc + 1) * step
        if args.eval_mode == "ball":
            z_center, z_radius = ellipse_arc_disk(theta_left, theta_right, rho)
            r = CBall(mid + half * z_center, abs(half) * z_radius)
            mat = contribution_ball(args.kind, mp.pi, s, r, n)
        elif args.arc_mode == "theta":
            z = ellipse_arc_box(theta_left, theta_right, rho)
            r_box = mp.iv.mpf([mid, mid]) + mp.iv.mpf([half, half]) * z
            mat = contribution_any(args.kind, mp.pi, s, r_box, n)
        else:
            z = ellipse_arc_rect_box(theta_left, theta_right, rho)
            r_box = mp.iv.mpf([mid, mid]) + mp.iv.mpf([half, half]) * z
            mat = contribution_any(args.kind, mp.pi, s, r_box, n)
        entry, row, frob = matrix_abs_bound(mat)
        if entry > max_entry:
            worst_arc = {
                "thetaLeft": f(theta_left),
                "thetaRight": f(theta_right),
                "entryBound": f(entry),
            }
        max_entry = max(max_entry, entry)
        max_row = max(max_row, row)
        max_frob = max(max_frob, frob)
    return max_entry, max_row, max_frob, worst_arc


def center_taylor_amplification(jet_order: int, needed_q: int) -> mp.mpf:
    max_amp = mp.mpf("0")
    for p in range(needed_q + 1):
        for i in range(jet_order):
            for j in range(jet_order):
                amp = mp.fsum(
                    mp.binomial(i + a, i)
                    * mp.binomial(j + p - a, j)
                    for a in range(p + 1)
                )
                max_amp = max(max_amp, amp)
    return max_amp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace-tail-json", default="endpoint_confluent_trace_tail_certificate.json")
    parser.add_argument(
        "--coefficient-json",
        default="endpoint_coefficient_interval_enclosure_consequence_theorem.json",
    )
    parser.add_argument("--kind", choices=["kb"], default="kb")
    parser.add_argument("--quadrature-order", type=int, default=200)
    parser.add_argument("--needed-trace-q", type=int, default=8)
    parser.add_argument("--segments", default="")
    parser.add_argument("--arcs", type=int, default=64)
    parser.add_argument("--segment-step", default="0.25")
    parser.add_argument("--eval-mode", choices=["majorant", "ball", "interval"], default="majorant")
    parser.add_argument("--arc-mode", choices=["rect", "theta"], default="rect")
    parser.add_argument("--rho-cap", default="2.0")
    parser.add_argument("--max-imag", default="0.20")
    parser.add_argument("--iv-dps", type=int, default=50)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--eigen-amplification", default="1e60")
    parser.add_argument("--s0", default="0.2825")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--max-trace-q", type=int, default=16)
    parser.add_argument("--matrix-rmax", default="12")
    parser.add_argument("--json-out", default="endpoint_confluent_segment_bernstein_certificate.json")
    args = parser.parse_args()

    if args.arcs < 4:
        raise SystemExit("--arcs must be at least 4")

    mp.mp.dps = args.dps
    mp.iv.dps = args.iv_dps
    trace_tail = load_json(args.trace_tail_json)
    coeff = load_json(args.coefficient_json)
    rho_cap = mp.mpf(args.rho_cap)
    max_imag = mp.mpf(args.max_imag)
    nodes = collocation_nodes(args)
    if args.segments.strip():
        segs = parse_mpf_list(args.segments)
    else:
        segs = refined_segments(mp.mpf(args.matrix_rmax), mp.mpf(args.segment_step))

    print(
        f"Endpoint confluent segment Bernstein certificate "
        f"nodes={len(nodes)} segments={len(segs)-1} order={args.quadrature_order}",
        flush=True,
    )

    rows = []
    total_entry_error = mp.mpf("0")
    total_row_error = mp.mpf("0")
    total_frob_error = mp.mpf("0")
    max_mrho = mp.mpf("0")
    worst = None
    for seg_idx, (a, b) in enumerate(zip(segs[:-1], segs[1:])):
        length = b - a
        half = length / 2
        rho = rho_from_imag_limit(half, max_imag, rho_cap)
        seg_entry_mrho = mp.mpf("0")
        seg_row_mrho = mp.mpf("0")
        seg_frob_mrho = mp.mpf("0")
        seg_worst = None
        for node_idx, s in enumerate(nodes):
            entry, row, frob, arc = segment_node_mrho(args, s, a, b, rho)
            if entry > seg_entry_mrho:
                seg_worst = {
                    "s": f(s),
                    "nodeIndex": node_idx,
                    "arc": arc,
                }
            seg_entry_mrho = max(seg_entry_mrho, entry)
            seg_row_mrho = max(seg_row_mrho, row)
            seg_frob_mrho = max(seg_frob_mrho, frob)
        entry_error = segment_error_bound(length, rho, args.quadrature_order, seg_entry_mrho)
        row_error = segment_error_bound(length, rho, args.quadrature_order, seg_row_mrho)
        frob_error = segment_error_bound(length, rho, args.quadrature_order, seg_frob_mrho)
        total_entry_error += entry_error
        total_row_error += row_error
        total_frob_error += frob_error
        if seg_entry_mrho > max_mrho:
            max_mrho = seg_entry_mrho
            worst = {
                "segmentIndex": seg_idx,
                "segment": [f(a), f(b)],
                "rho": f(rho),
                "mrhoEntry": f(seg_entry_mrho),
                "entryErrorBound": f(entry_error),
                "worstNodeArc": seg_worst,
            }
        rows.append(
            {
                "segmentIndex": seg_idx,
                "segment": [f(a), f(b)],
                "rho": f(rho),
                "mrhoEntry": f(seg_entry_mrho),
                "mrhoRow": f(seg_row_mrho),
                "mrhoFrobenius": f(seg_frob_mrho),
                "entryErrorBound": f(entry_error),
                "entryErrorBoundText": mp.nstr(entry_error, 16),
                "rowErrorBound": f(row_error),
                "frobeniusErrorBound": f(frob_error),
                "worstNodeArc": seg_worst,
            }
        )
        print(
            f"  seg {seg_idx+1:02d}/{len(segs)-1} [{fmt(a,5)},{fmt(b,5)}] "
            f"rho={fmt(rho,6)} entry_err={mp.nstr(entry_error, 8)}",
            flush=True,
        )

    d = args.jet_order
    amp_center = center_taylor_amplification(args.jet_order, args.needed_trace_q)
    center_entry_error = amp_center * total_entry_error
    center_spectral_error = d * center_entry_error
    eigen_amp = mp.mpf(args.eigen_amplification)
    trace_error_bound = eigen_amp * center_spectral_error
    trace_radius = mp.mpf(str(trace_tail["traceDerivativeEntryRadius"]))
    trace_target = mp.mpf(str(trace_tail["traceRadiusTarget"]))
    coeff_radius = min(
        mp.mpf(str(coeff["scaledCompanionIntervalRadius"])),
        mp.mpf(str(coeff["boundaryRowIntervalRadius"])),
    )
    declared_amplified_refinement_closed = trace_error_bound < trace_radius
    declared_amplified_target_closed = trace_error_bound < trace_target
    declared_amplified_coefficient_closed = trace_error_bound < coeff_radius
    center_refinement_closed = center_spectral_error < trace_radius
    center_target_closed = center_spectral_error < trace_target
    center_coefficient_closed = center_spectral_error < coeff_radius

    data = {
        "theoremName": "endpoint confluent segment Bernstein quadrature certificate",
        "traceTailJson": args.trace_tail_json,
        "coefficientJson": args.coefficient_json,
        "kind": args.kind,
        "quadratureOrder": args.quadrature_order,
        "segmentCount": len(segs) - 1,
        "segmentStep": args.segment_step,
        "segments": [[f(a), f(b)] for a, b in zip(segs[:-1], segs[1:])],
        "collocationNodeCount": len(nodes),
        "neededTraceQ": args.needed_trace_q,
        "bigOrder": args.jet_order + args.needed_trace_q,
        "arcs": args.arcs,
        "evalMode": args.eval_mode,
        "arcMode": args.arc_mode,
        "rhoCap": f(rho_cap),
        "maxImag": f(max_imag),
        "ivDps": args.iv_dps,
        "totalEntryQuadratureError": f(total_entry_error),
        "totalEntryQuadratureErrorText": mp.nstr(total_entry_error, 16),
        "totalRowQuadratureError": f(total_row_error),
        "totalFrobeniusQuadratureError": f(total_frob_error),
        "centerTaylorAmplificationMax": f(amp_center),
        "centerTaylorEntryErrorBound": f(center_entry_error),
        "centerTaylorEntryErrorBoundText": mp.nstr(center_entry_error, 16),
        "centerTaylorSpectralErrorBound": f(center_spectral_error),
        "eigenAmplificationBudget": f(eigen_amp),
        "traceDerivativeErrorBudget": f(trace_error_bound),
        "traceDerivativeErrorBudgetText": mp.nstr(trace_error_bound, 16),
        "traceDerivativeRefinementRadius": f(trace_radius),
        "traceDerivativeRefinementRadiusText": trace_tail.get("traceDerivativeEntryRadiusText"),
        "traceRadiusTarget": f(trace_target),
        "coefficientRadiusScale": f(coeff_radius),
        "quadratureBelowTraceRefinementRadius": bool(center_refinement_closed),
        "quadratureBelowTraceTarget": bool(center_target_closed),
        "quadratureBelowCoefficientScale": bool(center_coefficient_closed),
        "centerTaylorSpectralBelowTraceTarget": bool(center_target_closed),
        "centerTaylorSpectralBelowCoefficientScale": bool(center_coefficient_closed),
        "declaredEigenAmplifiedBelowTraceRefinementRadius": bool(declared_amplified_refinement_closed),
        "declaredEigenAmplifiedBelowTraceTarget": bool(declared_amplified_target_closed),
        "declaredEigenAmplifiedBelowCoefficientScale": bool(declared_amplified_coefficient_closed),
        "segmentRows": rows,
        "worstSegment": worst,
        "segmentBernsteinQuadratureStatus": status(
            "segment Bernstein Gauss-Legendre quadrature bound",
            bool(center_target_closed),
            (
                "Closed when the Bernstein/Chebyshev tail bound for the "
                "explicit Gauss-Legendre segments, propagated through the "
                "center Taylor map, is below the working trace target."
            ),
        ),
        "strictQuadratureReplacementStatus": status(
            "geometric quadrature-refinement replacement",
            bool(center_coefficient_closed),
            (
                "Closed when the deterministic segment Bernstein bound is "
                "small enough, at the displayed quadrature order and segment "
                "grid, to fit below the downstream finite Krawczyk coefficient "
                "radius scale before the separate eigenrow-formalization step."
            ),
        ),
        "remainingEigenrowFormalizationStatus": status(
            "formal interval eigenrow Taylor propagation",
            False,
            (
                "Open as a standalone theorem.  The script uses a declared "
                "eigen-amplification budget and the observed spectral gap; a "
                "fully formal proof would propagate the interval matrix "
                "through the eigenvector Taylor recurrence."
            ),
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("Endpoint confluent segment Bernstein certificate")
    print(f"  total entry quadrature error: {mp.nstr(total_entry_error, 12)}")
    print(f"  center spectral error bound: {mp.nstr(center_spectral_error, 12)}")
    print(f"  trace target: {mp.nstr(trace_target, 12)}")
    print(f"  coefficient radius scale: {mp.nstr(coeff_radius, 12)}")
    print(f"  center below trace target: {center_target_closed}")
    print(f"  center below coefficient scale: {center_coefficient_closed}")
    print(f"  declared eigen-amplified budget: {mp.nstr(trace_error_bound, 12)}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
