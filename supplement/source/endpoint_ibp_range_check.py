#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("endpoint_ibp_range_check.py requires numpy; run with python") from exc


def k_endpoint(c, x):
    return 2.0 * c * (x ** 1.25) * (2.0 * c * x - 3.0) * math.exp(-c * x)


def a_base(c, x):
    return (x ** 1.5) * math.exp(-c * x)


def h_kernel(c, r, s):
    # h_r(S) = -d/dS a(rS), positive for r,S >= 1 and c > 3/2.
    x = r * s
    return r * (c * x - 1.5) * math.sqrt(x) * math.exp(-c * x)


def sample_f(s):
    return math.exp(-0.7 * (s - 1.0)) * (1.0 - 0.22 * (s - 1.0) + 0.03 * (s - 1.0) ** 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--Smax", type=float, default=12.0)
    parser.add_argument("--Rmax", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=260)
    args = parser.parse_args()

    s_pts, s_wts = quadrature(args.Smax - 1.0, args.order)
    s_pts = s_pts + 1.0
    r_pts, r_wts = quadrature(args.Rmax - 1.0, args.order)
    r_pts = r_pts + 1.0

    direct = 0.0
    ibp = 0.0
    direct_deriv = 0.0
    ibp_deriv = 0.0
    min_h = float("inf")
    for r, rw in zip(r_pts, r_wts):
        t_val = 0.0
        h_val = 0.0
        dt_val = 0.0
        dh_val = 0.0
        for s, sw in zip(s_pts, s_wts):
            f = sample_f(float(s))
            x = float(r * s)
            kernel = k_endpoint(args.c, x)
            a = (s ** -0.25) * f
            h = h_kernel(args.c, float(r), float(s))
            min_h = min(min_h, h)
            t_val += float(sw) * kernel * f / float(s)
            h_val += float(sw) * 4.0 * args.c * (float(r) ** -0.25) * h * a
            dt_val += float(sw) * math.log(x) * kernel * f / float(s)
            dh_val += (
                float(sw)
                * 4.0
                * args.c
                * (float(r) ** -0.25)
                * math.log(x)
                * h
                * a
            )
        direct += float(rw) * t_val * t_val / float(r)
        ibp += float(rw) * h_val * h_val / float(r)
        direct_deriv += float(rw) * t_val * dt_val / float(r)
        ibp_deriv += float(rw) * h_val * dh_val / float(r)

    print(
        f"endpoint IBP range check S,R=[1,{args.Smax:g}] order={args.order} c={args.c:g}"
    )
    print(f"  min h_r(S) on grid       ={min_h:.12e}")
    print(f"  ||T F||^2 direct          ={direct:.12e}")
    print(f"  ||T F||^2 ibp             ={ibp:.12e}")
    print(f"  norm defect               ={direct - ibp:.12e}")
    print(f"  <TF, d_p TF> direct       ={direct_deriv:.12e}")
    print(f"  <TF, d_p TF> ibp          ={ibp_deriv:.12e}")
    print(f"  derivative defect         ={direct_deriv - ibp_deriv:.12e}")


if __name__ == "__main__":
    main()
