#!/usr/bin/env python3
import argparse
import math

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "second_order_cosh_ibp.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


PI = math.pi


def mode_right_derivative(n: int):
    c = PI * n * n
    total = 0.0
    for coeff, lam in (
        (4.0 * PI * PI * n * n * n * n, 4.5),
        (-6.0 * PI * n * n, 2.5),
    ):
        term = coeff * math.exp(-c)
        total += term * (lam - 2.0 * c)
    return total


def weights_for(kind: str):
    if kind == "raw1":
        return {1: 1.0}
    if kind == "raw2":
        return {1: 1.0, 2: 1.0}
    if kind == "raw3":
        return {1: 1.0, 2: 1.0, 3: 1.0}
    if kind == "tilde3":
        alpha3 = -(mode_right_derivative(1) + mode_right_derivative(2)) / mode_right_derivative(3)
        return {1: 1.0, 2: 1.0, 3: alpha3}
    raise ValueError(kind)


def simpson_grid(rmax: float, intervals: int):
    if intervals % 2:
        raise ValueError("intervals must be even")
    h = rmax / intervals
    xs = [h * i for i in range(intervals + 1)]
    ws = []
    for i in range(intervals + 1):
        factor = 1 if i == 0 or i == intervals else 4 if i % 2 else 2
        ws.append(h * factor / 3.0)
    return xs, ws


def finite_modes(kind: str):
    return [(int(n), float(weight)) for n, weight in sorted(weights_for(kind).items())]


def base_values(n: int, v: float):
    c = PI * n * n
    ev = math.exp(v)
    b = math.exp(0.25 * v - c * ev)
    # score a=-B'/B and multiplier m=(4D_v^2-1/4)B/B.
    a = c * ev - 0.25
    ap = c * ev
    m = 4.0 * c * c * ev * ev - 6.0 * c * ev
    mp = 8.0 * c * c * ev * ev - 6.0 * c * ev
    return b, a, ap, m, mp


def phi_from_base(n: int, t: float):
    v = 2.0 * t
    b, _, _, m, _ = base_values(n, v)
    return m * b


def phi_direct(n: int, t: float):
    c = PI * n * n
    e2 = math.exp(2.0 * t)
    return (4.0 * c * c * e2 * e2 - 6.0 * c * e2) * math.exp(0.5 * t - c * e2)


def weight_values(s: float, t: float, omega: float):
    # Same-sign kernel in v=2t variables has U=(s+t)/2.
    u = 0.5 * (s + t)
    ou = omega * u
    w = u * math.cosh(ou)
    # (partial_s + partial_t) is exactly d/du on functions of U.
    wp = math.cosh(ou) + omega * u * math.sinh(ou)
    return w, wp


def coeff_matrix(s: float, t: float, omega: float, modes, which: str):
    w, wp = weight_values(s, t, omega)
    left = [base_values(n, s) for n, _ in modes]
    right = [base_values(n, t) for n, _ in modes]
    out = np.zeros((len(modes), len(modes)), dtype=float)
    for i, (_, _) in enumerate(modes):
        _, ai, api, mi, mpi = left[i]
        for j, (_, _) in enumerate(modes):
            _, aj, apj, mj, mpj = right[j]
            q = ai + aj
            numerator = w * mi * mj
            c = numerator / q
            if which == "C":
                out[i, j] = c
            elif which == "D":
                out[i, j] = (
                    (wp * mi * mj + w * (mpi * mj + mi * mpj)) / q
                    - numerator * (api + apj) / (q * q)
                )
            elif which == "V":
                out[i, j] = numerator
            else:
                raise ValueError(which)
    return out


def scalar_base_values(v: float, modes):
    ev = math.exp(v)
    logs = []
    for n, weight in modes:
        c = PI * n * n
        logs.append(math.log(weight) + 0.25 * v - c * ev)
    shift = max(logs)
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    for n, weight in modes:
        c = PI * n * n
        log_base = math.log(weight) + 0.25 * v - c * ev
        base = math.exp(log_base - shift)
        p = 0.25 - c * ev
        p1 = -c * ev
        p2 = -c * ev
        s0 += base
        s1 += base * p
        s2 += base * (p * p + p1)
        s3 += base * (p * p * p + 3.0 * p * p1 + p2)
    psi = 4.0 * s2 - 0.25 * s0
    psi1 = 4.0 * s3 - 0.25 * s1
    a = -s1 / s0
    ap = -(s2 * s0 - s1 * s1) / (s0 * s0)
    mu = psi / s0
    mup = (psi1 * s0 - psi * s1) / (s0 * s0)
    log_b0 = shift + math.log(s0)
    b0 = math.exp(log_b0) if log_b0 > -745.0 else 0.0
    return b0, a, ap, mu, mup


def scalar_coeff(s: float, t: float, omega: float, modes, which: str):
    w, wp = weight_values(s, t, omega)
    _, a_s, ap_s, mu_s, mup_s = scalar_base_values(s, modes)
    _, a_t, ap_t, mu_t, mup_t = scalar_base_values(t, modes)
    q = a_s + a_t
    numerator = w * mu_s * mu_t
    if which == "C":
        return numerator / q
    if which == "D":
        return (
            (wp * mu_s * mu_t + w * (mup_s * mu_t + mu_s * mup_t)) / q
            - numerator * (ap_s + ap_t) / (q * q)
        )
    if which == "V":
        return numerator
    raise ValueError(which)


def scalar_kernel_value(x: float, y: float, omega: float, modes, which: str, u: float = 0.0):
    s = u + 2.0 * x
    t = u + 2.0 * y
    coeff = coeff_matrix(s, t, omega, modes, which)
    left = np.array([weight * base_values(n, s)[0] for n, weight in modes], dtype=float)
    right = np.array([weight * base_values(n, t)[0] for n, weight in modes], dtype=float)
    return float(left @ coeff @ right)


def scalar_base_kernel_value(
    x: float, y: float, omega: float, modes, which: str, u: float = 0.0
):
    s = u + 2.0 * x
    t = u + 2.0 * y
    b_s = scalar_base_values(s, modes)[0]
    b_t = scalar_base_values(t, modes)[0]
    return b_s * scalar_coeff(s, t, omega, modes, which) * b_t


def kernel_matrix(xs, omega: float, modes, which: str, u: float = 0.0):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = scalar_kernel_value(float(x), float(y), omega, modes, which, u)
            mat[i, j] = mat[j, i] = value
    return mat


def scalar_base_kernel_matrix(xs, omega: float, modes, which: str, u: float = 0.0):
    n = len(xs)
    mat = np.zeros((n, n), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = scalar_base_kernel_value(float(x), float(y), omega, modes, which, u)
            mat[i, j] = mat[j, i] = value
    return mat


def integrated_matrix(xs, omega: float, modes, which: str, umax: float, intervals: int):
    us, ws = simpson_grid(umax, intervals)
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for u, w in zip(us, ws):
        mat += float(w) * kernel_matrix(xs, omega, modes, which, float(u))
    return mat


def scalar_base_integrated_matrix(
    xs, omega: float, modes, which: str, umax: float, intervals: int
):
    us, ws = simpson_grid(umax, intervals)
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for u, w in zip(us, ws):
        mat += float(w) * scalar_base_kernel_matrix(xs, omega, modes, which, float(u))
    return mat


def direct_integrated_matrix(xs, omega: float, modes, umax: float, intervals: int):
    us, ws = simpson_grid(umax, intervals)
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for u, quad_w in zip(us, ws):
        for i, x in enumerate(xs):
            s = float(u) + 2.0 * float(x)
            left_phi = sum(weight * phi_from_base(n, 0.5 * s) for n, weight in modes)
            for j, y in enumerate(xs[: i + 1]):
                t = float(u) + 2.0 * float(y)
                right_phi = sum(weight * phi_from_base(n, 0.5 * t) for n, weight in modes)
                w, _ = weight_values(s, t, omega)
                value = float(quad_w) * w * left_phi * right_phi
                mat[i, j] += value
                if i != j:
                    mat[j, i] += value
    return mat


def eig_summary(name: str, mat):
    evals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
    print(f"  {name}: min={evals[0]:.12e} max={evals[-1]:.12e}")
    return evals[0]


def verify_identity(modes):
    worst = -1.0
    worst_at = None
    for n, _ in modes:
        for t in np.linspace(0.0, 4.0, 81):
            a = phi_from_base(n, float(t))
            b = phi_direct(n, float(t))
            err = abs(a - b) / max(1e-300, abs(b))
            if err > worst:
                worst = err
                worst_at = (n, float(t), a, b)
    print(
        "identity max rel error="
        f"{worst:.3e} at n={worst_at[0]} t={worst_at[1]:.6g}"
    )


def scan_vector_base(args, modes, xs):
    c0 = kernel_matrix(xs, args.omega, modes, "C", 0.0)
    eig_summary("boundary C(r=0)", c0)

    min_d_layer = float("inf")
    min_d_at = None
    for u in np.linspace(0.0, args.layer_umax, args.layer_n):
        d_layer = kernel_matrix(xs, args.omega, modes, "D", float(u))
        low = np.linalg.eigvalsh((d_layer + d_layer.T) / 2.0)[0]
        if low < min_d_layer:
            min_d_layer = float(low)
            min_d_at = float(u)
    print(f"  worst D layer: min={min_d_layer:.12e} at u={min_d_at:.6g}")

    d_int = integrated_matrix(xs, args.omega, modes, "D", args.umax, args.intervals)
    eig_summary("integrated D", d_int)

    # Constants 1/8 cancel in the comparison. The identity is
    # integral V = C(r=0) + integral D.
    direct = direct_integrated_matrix(xs, args.omega, modes, args.umax, args.intervals)
    reconstructed = c0 + d_int
    diff = direct - reconstructed
    scale = max(1e-300, float(np.max(np.abs(direct))))
    print(f"  IBP check max_abs={np.max(np.abs(diff)):.12e} rel={np.max(np.abs(diff)) / scale:.12e}")
    eig_summary("direct integral", direct)
    eig_summary("C + integrated D", reconstructed)


def scan_scalar_base(args, modes, xs):
    c0 = scalar_base_kernel_matrix(xs, args.omega, modes, "C", 0.0)
    eig_summary("scalar boundary C(r=0)", c0)

    min_d_layer = float("inf")
    min_d_at = None
    for u in np.linspace(0.0, args.layer_umax, args.layer_n):
        d_layer = scalar_base_kernel_matrix(xs, args.omega, modes, "D", float(u))
        low = np.linalg.eigvalsh((d_layer + d_layer.T) / 2.0)[0]
        if low < min_d_layer:
            min_d_layer = float(low)
            min_d_at = float(u)
    print(f"  scalar worst D layer: min={min_d_layer:.12e} at u={min_d_at:.6g}")

    d_int = scalar_base_integrated_matrix(xs, args.omega, modes, "D", args.umax, args.intervals)
    eig_summary("scalar integrated D", d_int)
    direct = direct_integrated_matrix(xs, args.omega, modes, args.umax, args.intervals)
    reconstructed = c0 + d_int
    diff = direct - reconstructed
    scale = max(1e-300, float(np.max(np.abs(direct))))
    print(
        "  scalar IBP check "
        f"max_abs={np.max(np.abs(diff)):.12e} rel={np.max(np.abs(diff)) / scale:.12e}"
    )
    eig_summary("direct integral", direct)
    eig_summary("scalar C + integrated D", reconstructed)


def scan(args):
    modes = finite_modes(args.kind)
    verify_identity(modes)
    xs = np.linspace(args.xmin, args.xmax, args.n)
    print(
        f"basis={args.basis} kind={args.kind} omega={args.omega:g} "
        f"x=[{args.xmin:g},{args.xmax:g}] "
        f"n={args.n} u=[0,{args.umax:g}] intervals={args.intervals}"
    )
    if args.basis == "vector":
        scan_vector_base(args, modes, xs)
    else:
        scan_scalar_base(args, modes, xs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--basis", choices=("vector", "scalar"), default="scalar")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--umax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--layer-umax", type=float, default=5.0)
    parser.add_argument("--layer-n", type=int, default=41)
    args = parser.parse_args()
    scan(args)


if __name__ == "__main__":
    main()
