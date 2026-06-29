#!/usr/bin/env python3
import argparse
import math

from klm_test import phi


def log_phi(x: float) -> float:
    return math.log(phi(x))


def score(x: float, h: float) -> float:
    return -(log_phi(x + h) - log_phi(x - h)) / (2.0 * h)


def score_prime(x: float, h: float) -> float:
    return -(log_phi(x + h) - 2.0 * log_phi(x) + log_phi(x - h)) / (h * h)


def c_base(s: float, t: float, omega: float, step: float) -> float:
    n = s + t
    hsum = score(s, step) + score(t, step)
    return 0.25 * n * math.cosh(omega * n) / hsum


def d_base(s: float, t: float, omega: float, step: float) -> float:
    n = s + t
    hsum = score(s, step) + score(t, step)
    hsum_prime = score_prime(s, step) + score_prime(t, step)
    w = 0.25 * n * math.cosh(omega * n)
    w_prime = 0.5 * math.cosh(omega * n) + 0.5 * omega * n * math.sinh(omega * n)
    return (w_prime * hsum - w * hsum_prime) / (hsum * hsum)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xmin", type=float, default=1e-3)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    xs = [args.xmin + (args.xmax - args.xmin) * i / (args.n - 1) for i in range(args.n)]
    min_c = float("inf")
    min_d = float("inf")
    max_d = -float("inf")
    min_d_at = None
    for s in xs:
        for t in xs:
            c_value = c_base(s, t, args.omega, args.h)
            d_value = d_base(s, t, args.omega, args.h)
            min_c = min(min_c, c_value)
            if d_value < min_d:
                min_d = d_value
                min_d_at = (s, t)
            max_d = max(max_d, d_value)

    print(f"omega={args.omega:g} grid=[{args.xmin:g},{args.xmax:g}] n={args.n}")
    print(f"  min C_base={min_c:.12e}")
    print(f"  min D_base={min_d:.12e} at s={min_d_at[0]:.6g}, t={min_d_at[1]:.6g}")
    print(f"  max D_base={max_d:.12e}")


if __name__ == "__main__":
    main()
