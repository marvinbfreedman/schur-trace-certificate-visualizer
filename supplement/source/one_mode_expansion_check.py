#!/usr/bin/env python3
import argparse
import math

from positive_branch_perturbation import quadrature

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("one_mode_expansion_check.py requires numpy; run with python") from exc


def b_value(p, c, x):
    return math.exp(p * x - c * math.exp(x))


def b_derivative(p, c, x):
    return (p - c * math.exp(x)) * b_value(p, c, x)


def hankel_apply(values, pts, wts, p, c, derivative=False):
    out = np.zeros_like(values)
    for i, r in enumerate(pts):
        total = 0.0
        for s, w, value in zip(pts, wts, values):
            kernel = b_derivative(p, c, float(r + s)) if derivative else b_value(
                p, c, float(r + s)
            )
            total += float(w) * kernel * value
        out[i] = total
    return out


def inner(f, g, wts):
    return float(np.sum(wts * f * g))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=1.495)
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--L", type=float, default=12.0)
    parser.add_argument("--order", type=int, default=260)
    args = parser.parse_args()

    pts, wts = quadrature(args.L, args.order)
    f = np.exp(-0.4 * pts) * (1.0 - 0.35 * pts + 0.05 * pts * pts)
    lf = pts * f
    a_const = 2.0 * args.p - 3.0

    y = hankel_apply(f, pts, wts, args.p, args.c)
    dy = hankel_apply(f, pts, wts, args.p, args.c, derivative=True)
    eta = hankel_apply(lf, pts, wts, args.p, args.c)
    deta = hankel_apply(lf, pts, wts, args.p, args.c, derivative=True)
    z = a_const * y - 2.0 * dy
    peta = a_const * eta - 2.0 * deta

    direct = inner(z, pts * z, wts) + inner(z, peta, wts)
    expanded = (
        a_const * a_const * (inner(y, pts * y, wts) + inner(y, eta, wts))
        + 2.0 * a_const * inner(y, y, wts)
        + 2.0 * a_const * y[0] * eta[0]
        + 4.0 * (inner(dy, pts * dy, wts) + inner(dy, deta, wts))
    )
    print(
        f"one-mode expansion check p={args.p:g} c={args.c:g} "
        f"L={args.L:g} order={args.order}"
    )
    print(f"  direct endpoint derivative form={direct:.12e}")
    print(f"  expanded formula            ={expanded:.12e}")
    print(f"  defect                      ={direct - expanded:.12e}")
    print("  expansion terms:")
    print(
        f"    A^2(<y,Ly>+<y,eta>)="
        f"{a_const * a_const * (inner(y, pts * y, wts) + inner(y, eta, wts)):.12e}"
    )
    print(f"    2A||y||^2                ={2.0 * a_const * inner(y, y, wts):.12e}")
    print(f"    2A y(0)eta(0)            ={2.0 * a_const * y[0] * eta[0]:.12e}")
    print(
        f"    4(<Dy,L Dy>+<Dy,D eta>)="
        f"{4.0 * (inner(dy, pts * dy, wts) + inner(dy, deta, wts)):.12e}"
    )


if __name__ == "__main__":
    main()
