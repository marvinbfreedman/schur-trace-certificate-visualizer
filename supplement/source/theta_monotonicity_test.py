#!/usr/bin/env python3
import argparse
import math

from klm_test import phi
from wedge_scalar_test import weight


def log_phi(x: float) -> float:
    value = phi(x)
    if value <= 0.0:
        return float("nan")
    return math.log(value)


def score(x: float, step: float) -> float:
    plus = log_phi(x + step)
    minus = log_phi(x - step)
    if not math.isfinite(plus) or not math.isfinite(minus):
        return float("nan")
    return -(plus - minus) / (2.0 * step)


def theta(c: float, v: float, r: float, omega: float, step: float) -> float:
    left = score(c - v + r, step)
    right = score(v + r, step)
    if not math.isfinite(left) or not math.isfinite(right):
        return float("nan")
    denom = left + right
    if denom == 0.0:
        return float("nan")
    return weight(c + 2.0 * r, omega) * (left - right) / denom


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--cmax", type=float, default=2.2)
    parser.add_argument("--rmax", type=float, default=0.4)
    parser.add_argument("--n-c", type=int, default=80)
    parser.add_argument("--n-r", type=int, default=80)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    min_boundary = float("inf")
    min_boundary_at = None
    min_step = float("inf")
    min_step_at = None

    for i in range(1, args.n_c):
        c = args.cmax * i / (args.n_c - 1)
        v_count = max(1, i // 2)
        for j in range(v_count + 1):
            v = 0.5 * c * j / v_count
            d = c - 2.0 * v
            psi0 = theta(c, v, 0.0, args.omega, args.h)
            if math.isfinite(psi0):
                boundary = psi0 - weight(d, args.omega)
                if boundary < min_boundary:
                    min_boundary = boundary
                    min_boundary_at = (c, v, psi0, weight(d, args.omega))

            previous = None
            previous_r = None
            for k in range(args.n_r):
                r = args.rmax * k / max(1, args.n_r - 1)
                value = theta(c, v, r, args.omega, args.h)
                if not math.isfinite(value):
                    continue
                if previous is not None:
                    diff = (value - previous) / (r - previous_r)
                    if diff < min_step:
                        min_step = diff
                        min_step_at = (c, v, previous_r, r, previous, value)
                previous = value
                previous_r = r

    print(f"omega={args.omega:g} cmax={args.cmax:g} rmax={args.rmax:g}")
    print(
        "  min boundary theta-W(d)="
        f"{min_boundary:.12e} at c={min_boundary_at[0]:.6g}, v={min_boundary_at[1]:.6g} "
        f"theta0={min_boundary_at[2]:.12e} Wd={min_boundary_at[3]:.12e}"
    )
    print(
        "  min theta r-slope="
        f"{min_step:.12e} at c={min_step_at[0]:.6g}, v={min_step_at[1]:.6g}, "
        f"r=[{min_step_at[2]:.6g},{min_step_at[3]:.6g}]"
    )


if __name__ == "__main__":
    main()
