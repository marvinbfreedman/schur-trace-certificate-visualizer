#!/usr/bin/env python3
r"""Between-node analytic enclosure of the endpoint negative eigenrow.

At a real center ``s0`` let ``J(s)`` denote the exact 9 by 9 confluent
endpoint matrix and let ``(e(s), lambda(s))`` be its normalized negative
eigenpair.  Nodewise interval eigenpair tests do not by themselves control
the open intervals between nodes.  This script supplies that missing local
argument.

The candidate eigenpair is its degree-N Taylor polynomial.  The exact Taylor
coefficients of J through degree N are enclosed by Arb quadrature of the
larger confluent matrix.  On a complex disk, a polynomial inverse of the
bordered eigenpair Jacobian is used in the Banach algebra H-infinity.  The
radii polynomial is

    Y + Z rho + 9 ||R|| rho^2 < rho.

The remaining Taylor tail of J is bounded by Cauchy's estimate on a larger
complex disk.  Once the radii polynomial closes, Cauchy's derivative estimate
gives a rigorous enclosure of e^(q), 0 <= q <= 8, on a strictly smaller real
cell.  A collection of such disks proves a between-node enclosure only when
their reported real cells cover the requested interval.

This is intentionally independent of the finite DAE collocation certificate.
It certifies the exact confluent integral/eigenpair field, not a fitted field.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import sys
from decimal import Decimal, localcontext
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402

from endpoint_uniform_eigenrow_continuation_theorem import (  # noqa: E402
    arb_lower_decimal,
    arb_upper_decimal,
    arb_upper_text,
    ball,
    clenshaw_curtis_rule,
    complex_s_disk_majorant,
    integrate_matrix_at_s,
    point,
    quadrature_error_bounds,
    real_tail_bound,
    series_div,
    series_exp,
    series_mul,
)
from global_trace_observability_gap import f  # noqa: E402
from trace_concomitant_exact_derivatives import eigen_taylor  # noqa: E402


DIM = 9
SYSTEM_DIM = DIM + 1


def status(label: str, closed: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": bool(closed),
        "status": "closed" if closed else "open",
        "reason": reason,
    }


def arb_abs_upper(value: arb) -> mp.mpf:
    return mp.mpf(str(arb_upper_decimal(abs(value))))


def arb_mid_mp(value: arb) -> mp.mpf:
    mid, _rad, exponent = value.mid().mid_rad_10exp()
    with localcontext() as context:
        context.prec = 110
        out = Decimal(int(mid)) * Decimal(10) ** int(exponent)
    return mp.mpf(str(out))


def arb_radius_upper(value: arb) -> mp.mpf:
    lo = mp.mpf(str(arb_lower_decimal(value)))
    hi = mp.mpf(str(arb_upper_decimal(value)))
    mid = arb_mid_mp(value)
    return max(mid - lo, hi - mid)


def mp_matrix_inf(matrix: mp.matrix) -> mp.mpf:
    return max(
        mp.fsum(abs(matrix[i, j]) for j in range(matrix.cols))
        for i in range(matrix.rows)
    )


def arb_matrix_inf_upper(matrix: list[list[arb]]) -> mp.mpf:
    return max(
        mp.fsum(arb_abs_upper(value) for value in row)
        for row in matrix
    )


def arb_vector_inf_upper(vector: list[arb]) -> mp.mpf:
    return max(arb_abs_upper(value) for value in vector)


def matrix_zero(rows: int, cols: int) -> list[list[arb]]:
    return [[arb(0) for _ in range(cols)] for _ in range(rows)]


def stable_endpoint_series(s: arb, r: arb, n: int) -> tuple[list[arb], list[arb]]:
    """Cancellation-preserving endpoint Taylor rows.

    The original quotient has denominator c(exp(s)-1).  Both quotients below
    are analytic at s=0; retaining ``exprel`` prevents the major overestimate
    caused by dividing numerator and denominator bounds separately.
    """
    c = arb.pi()
    x = r.exp()
    tau = x - 1
    one = [arb(0) for _ in range(n)]
    one[0] = 1
    s_series = [arb(0) for _ in range(n)]
    s_series[0] = s
    if n > 1:
        s_series[1] = 1
    exp_s = series_exp(s_series)
    delta = [c * value for value in exp_s]
    delta[0] -= c
    minus_delta_tau = [-tau * value for value in delta]
    exp_delta = series_exp(minus_delta_tau)
    exprel = series_exprel(minus_delta_tau)
    exp_c = (-c * tau).exp()

    vb = [
        exp_c
        * (
            -(1 - c * x) * tau * exprel[k]
            - x * exp_delta[k]
        )
        for k in range(n)
    ]

    lam = [c * value for value in exp_s]
    prefactor = list(one)
    for k in range(n):
        prefactor[k] -= x * lam[k]
    gb = [exp_c * value for value in series_mul(prefactor, exp_delta)]
    inv_exprel = series_div(one, series_exprel(s_series))
    wb = [value / c for value in series_mul(inv_exprel, gb)]
    return vb, wb


def series_exprel(values):
    """Taylor series for expm1(w)/w, valid even when w[0] contains zero."""
    n = len(values)
    if all(value == 0 for value in values):
        out = [arb(0) for _ in range(n)]
        out[0] = 1
        return out
    if isinstance(values[0], arb) and not values[0].contains(0):
        numerator = series_exp(values)
        numerator[0] -= 1
        return series_div(numerator, values)
    zero = acb(0) if isinstance(values[0], acb) else arb(0)
    one_scalar = acb(1) if isinstance(values[0], acb) else arb(1)
    u = list(values)
    base = values[0]
    u[0] = zero
    power = [zero for _ in range(n)]
    power[0] = one_scalar
    out = [zero for _ in range(n)]
    for k in range(n):
        derivative_over_factorial = (
            base.hypgeom_1f1(k + 1, k + 2)
            / point(k + 1)
            / math.factorial(k)
        )
        for j in range(n):
            out[j] += derivative_over_factorial * power[j]
        if k + 1 < n:
            power = series_mul(power, u)
    return out


def stable_contribution_matrix(s: arb, r: arb, n: int) -> list[list[arb]]:
    vb, wb = stable_endpoint_series(s, r, n)
    weight = (point("2.5") * r).exp()
    out = matrix_zero(n, n)
    for i in range(n):
        for j in range(i + 1):
            value = weight * (
                r * vb[i] * vb[j]
                + point("0.5") * (vb[i] * wb[j] + wb[i] * vb[j])
            )
            out[i][j] = value
            out[j][i] = value
    return out


def stable_segment_entry_majorant(
    center: Decimal,
    parameter_radius: Decimal,
    coefficient_radius: Decimal,
    a: Decimal,
    b: Decimal,
    rho: Decimal,
    n: int,
) -> arb:
    """Bernstein-ellipse majorant preserving the s=0 cancellation."""
    c = arb.pi()
    s0 = point(center)
    prad = point(parameter_radius)
    zrad = point(coefficient_radius)
    s_lo = s0 - prad
    s_hi = s0 + prad
    if arb_lower_decimal(s_lo) <= 0:
        raise ValueError("stable endpoint majorant requires center-parameter_radius>0")

    major = (point(rho) + 1 / point(rho)) / 2
    minor = (point(rho) - 1 / point(rho)) / 2
    mid = point((a + b) / 2)
    half = point((b - a) / 2)
    re_r_min = mid - half * major
    re_r_max = mid + half * major
    im_r = half * minor
    x_abs = re_r_max.exp()
    tau_abs = x_abs + 1
    re_tau_lo = re_r_min.exp() * im_r.cos() - 1
    im_tau_abs = re_r_max.exp() * abs(im_r.sin())

    exp_s_lo = s_lo.exp() * prad.cos()
    exp_s_hi = s_hi.exp()
    delta_re_lo = c * (exp_s_lo - 1)
    delta_re_hi = c * (exp_s_hi + 1)
    delta_im_abs = c * exp_s_hi * abs(prad.sin())
    delta_abs = c * (exp_s_hi + 1)
    re_tau_lower_decimal = arb_lower_decimal(re_tau_lo)
    delta_re_factor = delta_re_lo if re_tau_lower_decimal >= 0 else delta_re_hi
    re_delta_tau_lo = (
        delta_re_factor * re_tau_lo - delta_im_abs * im_tau_abs
    )
    re_lambda_tau_lo = c * re_tau_lo + re_delta_tau_lo

    exprel_bound = (-min(arb_lower_decimal(re_delta_tau_lo), Decimal(0))).exp()
    exp_delta_bound = (-arb_lower_decimal(re_delta_tau_lo)).exp()
    exp_c_bound = (-c * re_tau_lo).exp()
    vb_sup = exp_c_bound * (
        (1 + c * x_abs) * tau_abs * point(exprel_bound)
        + x_abs * point(exp_delta_bound)
    )

    endpoint_abs = abs(s0) + prad
    denominator_lower = s_lo.exp() - 1
    inv_exprel_bound = endpoint_abs / denominator_lower
    lam_abs = c * exp_s_hi
    gb_sup = (1 + x_abs * lam_abs) * (-re_lambda_tau_lo).exp()
    wb_sup = inv_exprel_bound * gb_sup / c

    r_abs = abs(mid) + abs(half) * major
    weight = (point("2.5") * re_r_max).exp()
    maximum = arb(0)
    powers = [arb(1)]
    for _ in range(1, n):
        powers.append(powers[-1] * zrad)
    for i in range(n):
        vi = vb_sup / powers[i]
        wi = wb_sup / powers[i]
        for j in range(i + 1):
            vj = vb_sup / powers[j]
            wj = wb_sup / powers[j]
            value = weight * (
                r_abs * vi * vj + point("0.5") * (vi * wj + wi * vj)
            )
            if arb_upper_decimal(value) > arb_upper_decimal(maximum):
                maximum = value
    return maximum


def stable_tail_bound(
    center: Decimal,
    parameter_radius: Decimal,
    coefficient_radius: Decimal,
    start_r: Decimal,
    n: int,
) -> arb:
    c = arb.pi()
    s0 = point(center)
    prad = point(parameter_radius)
    zrad = point(coefficient_radius)
    s_lo = s0 - prad
    s_hi = s0 + prad
    exp_s_lo = s_lo.exp() * prad.cos()
    exp_s_hi = s_hi.exp()
    lambda_re = c * exp_s_lo
    lambda_abs = c * exp_s_hi
    endpoint_abs = abs(s0) + prad
    inv_exprel = endpoint_abs / (s_lo.exp() - 1)
    cv = c.exp() * (point(2) + c)
    cw = inv_exprel * lambda_re.exp() * (1 + lambda_abs) / c
    coefficient_factor = zrad ** (-2 * (n - 1))
    constant = coefficient_factor * (cv * cv + cv * cw)
    x0 = point(start_r).exp()
    y = c * x0
    gamma_integer = math.factorial(7) * (-y).exp() * sum(
        (y**k / math.factorial(k) for k in range(8)), arb(0)
    )
    return constant * gamma_integer / c**8


def stable_quadrature_error_bounds(
    center: Decimal,
    parameter_radius: Decimal,
    coefficient_radius: Decimal,
    boundaries: list[Decimal],
    order: int,
    rho: Decimal,
    n: int,
) -> arb:
    total = arb(0)
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        major = (point(rho) + 1 / point(rho)) / 2
        minor = (point(rho) - 1 / point(rho)) / 2
        midpoint = point((a + b) / 2)
        half = point((b - a) / 2)
        r_box = acb(
            midpoint + ball(0, half * major),
            ball(0, half * minor),
        )
        matrix = stable_contribution_matrix(point(center), r_box, n)
        absolute_entries = [
            abs(matrix[i][j]) for i in range(n) for j in range(i + 1)
        ]
        if any(str(value) == "nan" for value in absolute_entries):
            majorant = stable_segment_entry_majorant(
                center,
                parameter_radius,
                coefficient_radius,
                a,
                b,
                rho,
                n,
            )
        else:
            majorant = max(absolute_entries, key=arb_upper_decimal)
        total += (
            point(8)
            * point(b - a)
            * majorant
            * point(rho) ** (-order)
            / (point(rho) - 1)
        )
    return total


def stable_integrate_matrix_at_s(
    s: arb,
    boundaries: list[Decimal],
    cc_nodes: list[arb],
    cc_weights: list[arb],
    n: int,
    common_error: arb,
) -> list[list[arb]]:
    out = matrix_zero(n, n)
    for a_dec, b_dec in zip(boundaries[:-1], boundaries[1:]):
        mid = point((a_dec + b_dec) / 2)
        half = point((b_dec - a_dec) / 2)
        for x, weight in zip(cc_nodes, cc_weights):
            r = mid + half * x
            matrix = stable_contribution_matrix(s, r, n)
            factor = half * weight
            for i in range(n):
                for j in range(i + 1):
                    out[i][j] += factor * matrix[i][j]
    error_ball = ball(0, arb_upper_text(common_error))
    for i in range(n):
        for j in range(i + 1):
            out[i][j] += error_ball
            out[j][i] = out[i][j]
    return out


def stable_complex_s_disk_majorant(
    center: Decimal,
    disk_radius: Decimal,
    step: Decimal,
    n: int,
    r_cutoff: Decimal = Decimal(6),
) -> arb:
    clearance = center - disk_radius
    if clearance <= 0:
        raise ValueError("complex disk crosses the endpoint center s=0")
    s_box = acb(ball(center, disk_radius), ball(0, disk_radius))
    auxiliary = min(Decimal("0.01"), clearance / 2)
    total = arb(0)
    a = Decimal(0)
    while a < r_cutoff:
        b = min(r_cutoff, a + step)
        r_box = ball((a + b) / 2, (b - a) / 2)
        matrix = stable_contribution_matrix(s_box, r_box, n)
        absolute_entries = [
            abs(matrix[i][j]) for i in range(n) for j in range(i + 1)
        ]
        if any(str(value) == "nan" for value in absolute_entries):
            majorant = stable_segment_entry_majorant(
                center,
                disk_radius + auxiliary,
                auxiliary,
                a,
                b,
                Decimal(1),
                n,
            )
        else:
            majorant = max(absolute_entries, key=arb_upper_decimal)
        total += point(b - a) * majorant
        a = b
    return total + stable_tail_bound(
        center, disk_radius + auxiliary, auxiliary, r_cutoff, n
    )


def matrix_series_mul(
    left: list[mp.matrix], right: list[list[list[arb]]]
) -> list[list[list[arb]]]:
    degree = len(left) + len(right) - 2
    out = [matrix_zero(left[0].rows, right[0][0].__len__()) for _ in range(degree + 1)]
    rows = left[0].rows
    inner = left[0].cols
    cols = len(right[0][0])
    for n in range(degree + 1):
        for p in range(max(0, n - len(right) + 1), min(n, len(left) - 1) + 1):
            q = n - p
            for i in range(rows):
                for j in range(cols):
                    out[n][i][j] += sum(
                        (point(left[p][i, k]) * right[q][k][j] for k in range(inner)),
                        arb(0),
                    )
    return out


def series_sup_matrix(coefficients: list[list[list[arb]]], radius: mp.mpf) -> mp.mpf:
    return mp.fsum(
        arb_matrix_inf_upper(coefficient) * radius**n
        for n, coefficient in enumerate(coefficients)
    )


def series_sup_vector(coefficients: list[list[arb]], radius: mp.mpf) -> mp.mpf:
    return mp.fsum(
        arb_vector_inf_upper(coefficient) * radius**n
        for n, coefficient in enumerate(coefficients)
    )


def load_candidate(cache_path: str, degree: int) -> tuple[mp.mpf, list[list[mp.mpf]], list[mp.mpf]]:
    data = json.loads(Path(cache_path).read_text(encoding="utf-8"))
    e_derivs = data["eDerivs"]
    lambda_derivs = data["lambdaDerivs"]
    if len(e_derivs) <= degree or len(lambda_derivs) <= degree:
        raise ValueError(f"{cache_path} does not contain a degree-{degree} jet")
    e = [
        [mp.mpf(str(value)) / mp.factorial(q) for value in e_derivs[q]]
        for q in range(degree + 1)
    ]
    lam = [
        mp.mpf(str(lambda_derivs[q])) / mp.factorial(q)
        for q in range(degree + 1)
    ]
    return mp.mpf(str(data.get("s", data.get("key", {}).get("s")))), e, lam


def load_cache_center(cache_path: str) -> mp.mpf:
    data = json.loads(Path(cache_path).read_text(encoding="utf-8"))
    return mp.mpf(str(data.get("s", data.get("key", {}).get("s"))))


def discover_candidates(cache_dir: str, degree: int = 0) -> list[tuple[mp.mpf, str]]:
    rows: dict[str, tuple[mp.mpf, str]] = {}
    for path in glob.glob(str(Path(cache_dir) / "*.json")):
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
            if len(data.get("eDerivs", [])) <= degree:
                continue
            s_text = str(data.get("s", data.get("key", {}).get("s")))
            s = mp.mpf(s_text)
            key = mp.nstr(s, 40)
            rows.setdefault(key, (s, path))
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
    return sorted(rows.values(), key=lambda row: row[0])


def exact_taylor_matrix_balls(
    center: mp.mpf,
    degree: int,
    quadrature_order: int,
    r_cutoff: Decimal,
    r_step: Decimal,
    cache_dir: str | None = None,
) -> tuple[list[list[list[arb]]], dict]:
    big_n = DIM + degree
    cache_path = None
    cache_key = {
        "center": mp.nstr(center, 80),
        "degree": degree,
        "quadratureOrder": quadrature_order,
        "rCutoff": str(r_cutoff),
        "rStep": str(r_step),
        "arbBits": int(ctx.prec),
    }
    if center < mp.mpf("0.08"):
        cache_key["formulaVersion"] = "stable-exprel-v2"
    if cache_dir:
        directory = Path(cache_dir)
        directory.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256(json.dumps(cache_key, sort_keys=True).encode()).hexdigest()[:24]
        cache_path = directory / f"exact_taylor_{digest}.json"
        if cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            coefficients = [
                [[arb(value) for value in row] for row in matrix]
                for matrix in cached["coefficients"]
            ]
            metadata = dict(cached["metadata"])
            metadata["cacheHit"] = True
            metadata["cachePath"] = str(cache_path)
            return coefficients, metadata
    boundaries = [Decimal(0)]
    while boundaries[-1] < r_cutoff:
        boundaries.append(min(r_cutoff, boundaries[-1] + r_step))
    center_dec = Decimal(mp.nstr(center, 80))
    stable_radius = min(Decimal("0.01"), center_dec * Decimal("0.45"))
    use_stable = center_dec < Decimal("0.08")
    if use_stable:
        quadrature_error = stable_quadrature_error_bounds(
            center_dec,
            stable_radius,
            stable_radius,
            boundaries,
            quadrature_order,
            Decimal(2),
            big_n,
        )
        tail_error = stable_tail_bound(
            center_dec, stable_radius, stable_radius, r_cutoff, big_n
        )
    else:
        quadrature_error, _rows = quadrature_error_bounds(
            center_dec,
            center_dec,
            boundaries,
            quadrature_order,
            Decimal(2),
            big_n,
        )
        tail_error = real_tail_bound(center_dec, center_dec, r_cutoff, big_n)
    common_error = quadrature_error + tail_error
    nodes, weights = clenshaw_curtis_rule(quadrature_order)
    integration_function = stable_integrate_matrix_at_s if use_stable else integrate_matrix_at_s
    big = integration_function(point(center_dec), boundaries, nodes, weights, big_n, common_error)
    coefficients = []
    for p in range(degree + 1):
        matrix = matrix_zero(DIM, DIM)
        for i in range(DIM):
            for j in range(DIM):
                matrix[i][j] = sum(
                    (
                        point(math.comb(i + a, i))
                        * point(math.comb(j + p - a, j))
                        * big[i + a][j + p - a]
                        for a in range(p + 1)
                    ),
                    arb(0),
                )
        coefficients.append(matrix)
    metadata = {
        "quadratureOrder": quadrature_order,
        "rCutoff": str(r_cutoff),
        "rStep": str(r_step),
        "quadratureEntryErrorUpper": arb_upper_text(quadrature_error),
        "rTailEntryErrorUpper": arb_upper_text(tail_error),
        "commonBigMatrixEntryErrorUpper": arb_upper_text(common_error),
        "cancellationPreservingEndpointFormula": use_stable,
        "stableCauchyRadius": str(stable_radius) if use_stable else None,
        "maxTaylorCoefficientEntryRadius": f(
            max(
                arb_radius_upper(value)
                for matrix in coefficients
                for row in matrix
                for value in row
            )
        ),
        "cacheHit": False,
        "cachePath": str(cache_path) if cache_path is not None else None,
    }
    if cache_path is not None:
        cache_path.write_text(
            json.dumps(
                {
                    "key": cache_key,
                    "coefficients": [
                        [[str(value) for value in row] for row in matrix]
                        for matrix in coefficients
                    ],
                    "metadata": metadata,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    return coefficients, metadata


def candidate_jacobian_series(
    a_coefficients: list[list[list[arb]]],
    e: list[list[mp.mpf]],
    lam: list[mp.mpf],
) -> list[list[list[arb]]]:
    degree = len(e) - 1
    out = [matrix_zero(SYSTEM_DIM, SYSTEM_DIM) for _ in range(degree + 1)]
    for n in range(degree + 1):
        for i in range(DIM):
            for j in range(DIM):
                out[n][i][j] = a_coefficients[n][i][j]
                if i == j:
                    out[n][i][j] -= point(lam[n])
            out[n][i][DIM] = -point(e[n][i])
            out[n][DIM][i] = point(e[n][i])
    return out


def candidate_residual_series(
    a_coefficients: list[list[list[arb]]],
    e: list[list[mp.mpf]],
    lam: list[mp.mpf],
) -> list[list[arb]]:
    degree = len(e) - 1
    out = [[arb(0) for _ in range(SYSTEM_DIM)] for _ in range(2 * degree + 1)]
    for n in range(2 * degree + 1):
        for p in range(max(0, n - degree), min(n, degree) + 1):
            q = n - p
            if p <= degree:
                for i in range(DIM):
                    out[n][i] += sum(
                        (a_coefficients[p][i][j] * point(e[q][j]) for j in range(DIM)),
                        arb(0),
                    )
            for i in range(DIM):
                out[n][i] -= point(lam[p]) * point(e[q][i])
            out[n][DIM] += point("0.5") * sum(
                (point(e[p][j]) * point(e[q][j]) for j in range(DIM)), arb(0)
            )
        if n == 0:
            out[n][DIM] -= point("0.5")
    return out


def midpoint_matrix(matrix: list[list[arb]]) -> mp.matrix:
    return mp.matrix([[arb_mid_mp(value) for value in row] for row in matrix])


def candidate_from_matrix_balls(
    coefficients: list[list[list[arb]]], degree: int
) -> tuple[list[list[mp.mpf]], list[mp.mpf]]:
    midpoint_coefficients = [midpoint_matrix(matrix) for matrix in coefficients]
    _values, e_derivs, lambda_derivs = eigen_taylor(midpoint_coefficients, degree)
    e = [
        [value / mp.factorial(q) for value in e_derivs[q]]
        for q in range(degree + 1)
    ]
    lam = [lambda_derivs[q] / mp.factorial(q) for q in range(degree + 1)]
    return e, lam


def inverse_series(jacobian: list[list[list[arb]]]) -> list[mp.matrix]:
    degree = len(jacobian) - 1
    h = [midpoint_matrix(matrix) for matrix in jacobian]
    out = [h[0] ** -1]
    for n in range(1, degree + 1):
        total = mp.matrix(SYSTEM_DIM)
        for k in range(1, n + 1):
            total += h[k] * out[n - k]
        out.append(-out[0] * total)
    return out


def identity_defect_series(
    inverse: list[mp.matrix], jacobian: list[list[list[arb]]]
) -> list[list[list[arb]]]:
    product = matrix_series_mul(inverse, jacobian)
    for i in range(SYSTEM_DIM):
        product[0][i][i] -= 1
    return product


def cauchy_tail(majorant: mp.mpf, inner: mp.mpf, outer: mp.mpf, degree: int) -> mp.mpf:
    ratio = inner / outer
    return 4 * majorant * ratio ** (degree + 1) / (1 - ratio)


def polynomial_derivative_interval(
    coefficients: list[mp.mpf], q: int, left: mp.mpf, right: mp.mpf
) -> arb:
    x = ball((left + right) / 2, (right - left) / 2)
    powers = [arb(1)]
    for _ in range(1, len(coefficients) - q):
        powers.append(powers[-1] * x)
    return sum(
        (
            point(math.factorial(n) // math.factorial(n - q) * coefficients[n])
            * powers[n - q]
            for n in range(q, len(coefficients))
        ),
        arb(0),
    )


def certify_disk(
    center: mp.mpf,
    cache_path: str | None,
    degree: int,
    inner_radius: mp.mpf,
    real_half_width: mp.mpf,
    outer_radius: mp.mpf,
    quadrature_order: int,
    r_cutoff: Decimal,
    r_step: Decimal,
    majorant_step: Decimal,
    matrix_cache_dir: str | None,
) -> dict:
    if not (0 < real_half_width < inner_radius < outer_radius < center):
        raise ValueError(
            f"disk radii at {center} must satisfy 0 < cell < inner < outer < center"
        )
    a_coefficients, quadrature = exact_taylor_matrix_balls(
        center, degree, quadrature_order, r_cutoff, r_step, matrix_cache_dir
    )
    candidate_source = "exact Arb Taylor-matrix midpoints"
    if cache_path is not None:
        try:
            cached_center, e, lam = load_candidate(cache_path, degree)
            if abs(cached_center - center) > mp.mpf("1e-40"):
                raise ValueError("cache center mismatch")
            candidate_source = cache_path
        except ValueError:
            e, lam = candidate_from_matrix_balls(a_coefficients, degree)
    else:
        e, lam = candidate_from_matrix_balls(a_coefficients, degree)
    jacobian = candidate_jacobian_series(a_coefficients, e, lam)
    inverse = inverse_series(jacobian)
    defect = identity_defect_series(inverse, jacobian)
    residual = candidate_residual_series(a_coefficients, e, lam)

    center_dec = Decimal(mp.nstr(center, 80))
    outer_dec = Decimal(mp.nstr(outer_radius, 80))
    majorant_ball = (
        stable_complex_s_disk_majorant(center_dec, outer_dec, majorant_step, DIM)
        if center_dec < Decimal("0.08")
        else complex_s_disk_majorant(center_dec, outer_dec, majorant_step, DIM)
    )
    majorant = mp.mpf(str(arb_upper_decimal(majorant_ball)))
    matrix_tail = cauchy_tail(majorant, inner_radius, outer_radius, degree)

    inverse_sup = mp.fsum(
        mp_matrix_inf(matrix) * inner_radius**n for n, matrix in enumerate(inverse)
    )
    e_sup = max(
        mp.fsum(abs(e[n][i]) * inner_radius**n for n in range(degree + 1))
        for i in range(DIM)
    )
    residual_sup = series_sup_vector(residual, inner_radius) + DIM * matrix_tail * e_sup
    defect_sup = series_sup_matrix(defect, inner_radius) + inverse_sup * DIM * matrix_tail
    y_bound = inverse_sup * residual_sup
    z_bound = defect_sup
    quadratic = DIM * inverse_sup
    discriminant = (1 - z_bound) ** 2 - 4 * quadratic * y_bound
    radius = mp.inf
    if discriminant > 0 and z_bound < 1:
        small_root = 2 * y_bound / ((1 - z_bound) + mp.sqrt(discriminant))
        # The radii polynomial vanishes at ``small_root``.  Select a point
        # strictly inside its negative interval rather than testing the root
        # itself; the factor two is valid whenever it remains below the large
        # root, which is checked by the inequality below.
        radius = 2 * small_root
    contraction = z_bound + quadratic * radius if mp.isfinite(radius) else mp.inf
    radii_closed = (
        mp.isfinite(radius)
        and radius > 0
        and y_bound + z_bound * radius + quadratic * radius**2 < radius
        and contraction < 1
    )

    margin = inner_radius - real_half_width
    derivative_errors = [
        mp.factorial(q) * radius / margin**q if radii_closed else mp.inf
        for q in range(9)
    ]
    cell_left = center - real_half_width
    cell_right = center + real_half_width
    derivative_rows = []
    for q in range(9):
        component_bounds = []
        for i in range(DIM):
            candidate = polynomial_derivative_interval(
                [e[n][i] for n in range(degree + 1)],
                q,
                -real_half_width,
                real_half_width,
            )
            enclosed = candidate + ball(0, derivative_errors[q])
            lower = arb_lower_decimal(enclosed)
            upper = arb_upper_decimal(enclosed)
            component_bounds.append(
                {
                    "lower": f(mp.mpf(str(lower))),
                    "upper": f(mp.mpf(str(upper))),
                    "lowerText": format(lower, ".80E"),
                    "upperText": format(upper, ".80E"),
                }
            )
        derivative_rows.append(
            {
                "q": q,
                "cauchyErrorRadius": f(derivative_errors[q]),
                "components": component_bounds,
            }
        )
    return {
        "center": f(center),
        "candidateSource": candidate_source,
        "orientationCachePath": cache_path,
        "degree": degree,
        "innerDiskRadius": f(inner_radius),
        "outerDiskRadius": f(outer_radius),
        "realCell": [f(cell_left), f(cell_right)],
        "realCellHalfWidth": f(real_half_width),
        "cauchyMargin": f(margin),
        "quadrature": quadrature,
        "complexMatrixMajorantUpper": f(majorant),
        "matrixTaylorTailUpper": f(matrix_tail),
        "polynomialInverseSup": f(inverse_sup),
        "residualSup": f(residual_sup),
        "defectSup": f(defect_sup),
        "radiiY": f(y_bound),
        "radiiZ": f(z_bound),
        "radiiQuadraticCoefficient": f(quadratic),
        "radiiDiscriminant": f(discriminant),
        "solutionRadius": f(radius),
        "contractionBound": f(contraction),
        "radiiPolynomialClosed": bool(radii_closed),
        "derivativeEnclosures": derivative_rows,
    }


def interval_covered(rows: list[dict], left: mp.mpf, right: mp.mpf) -> tuple[bool, mp.mpf]:
    intervals = sorted(
        (mp.mpf(str(row["realCell"][0])), mp.mpf(str(row["realCell"][1])))
        for row in rows
        if row["radiiPolynomialClosed"]
    )
    frontier = left
    for a, b in intervals:
        if a > frontier:
            return False, frontier
        frontier = max(frontier, b)
        if frontier >= right:
            return True, frontier
    return False, frontier


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", default="0.02")
    parser.add_argument("--right", default="0.545")
    parser.add_argument("--cache-dir", default=".endpoint_flow_cache")
    parser.add_argument(
        "--matrix-cache-dir", default=".endpoint_between_node_cauchy_cache"
    )
    parser.add_argument("--cache", action="append", default=[])
    parser.add_argument("--degree", type=int, default=16)
    parser.add_argument("--center", action="append", default=[])
    parser.add_argument("--inner-radius", default="0.006")
    parser.add_argument("--real-half-width", default="0.0033")
    parser.add_argument("--outer-radius", default="0.3")
    parser.add_argument("--quadrature-order", type=int, default=192)
    parser.add_argument("--r-cutoff", default="6")
    parser.add_argument("--r-step", default="0.25")
    parser.add_argument("--majorant-step", default="0.02")
    parser.add_argument("--bits", type=int, default=256)
    parser.add_argument(
        "--json-out", default="endpoint_eigenrow_between_node_cauchy_certificate.json"
    )
    args = parser.parse_args()

    ctx.prec = args.bits
    mp.mp.dps = max(90, int(args.bits * 0.32))
    left = mp.mpf(args.left)
    right = mp.mpf(args.right)
    available = discover_candidates(args.cache_dir)
    work: list[tuple[mp.mpf, str | None]] = []
    for cache_path in args.cache:
        work.append((load_cache_center(cache_path), cache_path))
    if args.center:
        for text_center in args.center:
            target = mp.mpf(text_center)
            nearest = min(available, key=lambda row: abs(row[0] - target)) if available else None
            cache_path = (
                nearest[1]
                if nearest is not None and abs(nearest[0] - target) < mp.mpf("1e-40")
                else None
            )
            work.append((target, cache_path))
    if not work:
        # A one-disk default is deliberate: it makes an expensive first run a
        # local validation rather than silently claiming global coverage.
        if not available:
            raise SystemExit("no candidate cache found and no --center supplied")
        nearest = min(available, key=lambda row: abs(row[0] - right))
        work = [nearest]

    rows = []
    for index, (center, cache_path) in enumerate(work, 1):
        outer_fraction = mp.mpf("0.95") if center < mp.mpf("0.08") else mp.mpf("0.75")
        outer = min(mp.mpf(args.outer_radius), center * outer_fraction)
        inner = min(mp.mpf(args.inner_radius), outer * mp.mpf("0.35"))
        half = min(mp.mpf(args.real_half_width), inner * mp.mpf("0.55"))
        print(
            f"disk {index}/{len(work)} center={mp.nstr(center, 14)} "
            f"cell={mp.nstr(half, 6)} inner={mp.nstr(inner, 6)} outer={mp.nstr(outer, 6)}",
            flush=True,
        )
        row = certify_disk(
            center,
            cache_path,
            args.degree,
            inner,
            half,
            outer,
            args.quadrature_order,
            Decimal(args.r_cutoff),
            Decimal(args.r_step),
            Decimal(args.majorant_step),
            args.matrix_cache_dir,
        )
        rows.append(row)
        print(
            f"  closed={row['radiiPolynomialClosed']} "
            f"Y={row['radiiY']:.3e} Z={row['radiiZ']:.3e} "
            f"rho={row['solutionRadius']:.3e}",
            flush=True,
        )

    covered, frontier = interval_covered(rows, left, right)
    local_closed = all(row["radiiPolynomialClosed"] for row in rows)
    exact_closed = local_closed and covered
    data = {
        "schemaVersion": 1,
        "theoremName": "between-node Bernstein/Cauchy endpoint eigenrow enclosure",
        "proofClass": "interval/ball certificate",
        "statement": (
            "On every reported closed real cell, the exact normalized negative eigenrow "
            "of the confluent endpoint matrix lies in the displayed analytic tube and its "
            "derivatives e^(q), 0<=q<=8, lie in the displayed Cauchy enclosures. The full "
            "interval statement is asserted only if the closed cells cover [left,right]."
        ),
        "domain": {
            "parameterInterval": [f(left), f(right)],
            "matrixFamily": "exact 9x9 confluent endpoint integral J(s)",
            "normalization": "e(s)^T e(s)=1 with orientation inherited from the cached negative eigenrow",
            "functionSpace": "bounded analytic eigenpair functions on each certified complex disk",
        },
        "dependencies": [
            "endpoint_confluent_segment_bernstein_certificate.json",
            "endpoint_eigenrow_interval_propagation_200.json",
        ],
        "method": {
            "matrixCoefficients": "Arb Clenshaw-Curtis quadrature plus Bernstein segment error and analytic r-tail",
            "betweenNode": "H-infinity Newton-Krawczyk with polynomial inverse Jacobian",
            "matrixTail": "Cauchy estimate from a larger exact complex-s disk majorant",
            "derivatives": "Cauchy derivative estimate on the cell-to-disk margin",
        },
        "disks": rows,
        "coverageFrontier": f(frontier),
        "localDiskCertificatesClosed": local_closed,
        "requestedIntervalCovered": covered,
        "betweenNodeEigenrowDerivativeEnclosureClosed": exact_closed,
        "betweenNodeEigenrowDerivativeStatus": status(
            "exact between-node e^(q), q<=8 enclosure",
            exact_closed,
            (
                "Every analytic radii polynomial closes and the certified real cells cover "
                "the requested interval."
                if exact_closed
                else "Only the listed local cells are certified; full interval coverage is not yet established."
            ),
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"requested interval covered: {covered} (frontier={mp.nstr(frontier, 14)})")
    print(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
