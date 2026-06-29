#!/usr/bin/env python3
import argparse
import math

from second_order_cosh_ibp import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "reduced_absorption_raw1.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


C = math.pi


def vals(v):
    ev = math.exp(v)
    a = C * ev - 0.25
    ap = C * ev
    mu = 4.0 * C * C * ev * ev - 6.0 * C * ev
    mup = 8.0 * C * C * ev * ev - 6.0 * C * ev
    log_b = 0.25 * v - C * ev
    return a, ap, mu, mup, log_b


def weight(s, t, omega):
    u = 0.5 * (s + t)
    ou = omega * u
    w = u * math.cosh(ou)
    wp = math.cosh(ou) + omega * u * math.sinh(ou)
    return w, wp


def boundary_reduced(s, t, omega):
    a_s, _, _, _, _ = vals(s)
    a_t, _, _, _, _ = vals(t)
    w, _ = weight(s, t, omega)
    return w / (a_s + a_t)


def d_coeff(s, t, omega):
    a_s, ap_s, mu_s, mup_s, _ = vals(s)
    a_t, ap_t, mu_t, mup_t, _ = vals(t)
    w, wp = weight(s, t, omega)
    q = a_s + a_t
    numerator = w * mu_s * mu_t
    return (
        (wp * mu_s * mu_t + w * (mup_s * mu_t + mu_s * mup_t)) / q
        - numerator * (ap_s + ap_t) / (q * q)
    )


def tail_reduced(s, t, omega, umax, intervals):
    _, _, mu_s, _, log_b_s = vals(s)
    _, _, mu_t, _, log_b_t = vals(t)
    us, ws = simpson_grid(umax, intervals)
    total = 0.0
    for u, qw in zip(us, ws):
        su = s + u
        tu = t + u
        *_, log_b_su = vals(su)
        *_, log_b_tu = vals(tu)
        log_ratio = log_b_su - log_b_s + log_b_tu - log_b_t
        if log_ratio < -745.0:
            continue
        total += qw * d_coeff(su, tu, omega) * math.exp(log_ratio) / (mu_s * mu_t)
    return total


def matrix(points, omega, which, umax, intervals):
    n = len(points)
    out = np.zeros((n, n), dtype=float)
    for i, s in enumerate(points):
        for j, t in enumerate(points[: i + 1]):
            if which == "boundary":
                value = boundary_reduced(float(s), float(t), omega)
            elif which == "tail":
                value = tail_reduced(float(s), float(t), omega, umax, intervals)
            elif which == "full":
                value = boundary_reduced(float(s), float(t), omega) + tail_reduced(
                    float(s), float(t), omega, umax, intervals
                )
            else:
                raise ValueError(which)
            out[i, j] = out[j, i] = value
    return out


def points_for(args):
    if args.grid == "linear":
        return np.linspace(args.vmin, args.vmax, args.n)
    if args.grid == "quadratic":
        z = np.linspace(0.0, 1.0, args.n)
        return args.vmin + (args.vmax - args.vmin) * z * z
    if args.grid == "geometric":
        if args.vmin <= 0:
            raise ValueError("geometric grid requires --vmin > 0")
        return np.geomspace(args.vmin, args.vmax, args.n)
    raise ValueError(args.grid)


def eigs(mat):
    return np.linalg.eigvalsh((mat + mat.T) / 2.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=5.2)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--umax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=320)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    pts = points_for(args)
    boundary = matrix(pts, args.omega, "boundary", args.umax, args.intervals)
    tail = matrix(pts, args.omega, "tail", args.umax, args.intervals)
    full = boundary + tail
    bvals = eigs(boundary)
    tvals = eigs(tail)
    fvals = eigs(full)
    print(
        f"raw1 reduced absorption omega={args.omega:g} grid={args.grid} "
        f"v=[{args.vmin:g},{args.vmax:g}] n={args.n} u=[0,{args.umax:g}]"
    )
    for name, vals_ in (("boundary", bvals), ("tail", tvals), ("full", fvals)):
        print(
            f"  {name}: neg={(vals_ < -args.tol).sum()} "
            f"zero={(abs(vals_) <= args.tol).sum()} pos={(vals_ > args.tol).sum()} "
            f"min={vals_[0]:.12e} max={vals_[-1]:.12e}"
        )
        print("    low=", " ".join(f"{v:.4e}" for v in vals_[: min(8, len(vals_))]))


if __name__ == "__main__":
    main()
