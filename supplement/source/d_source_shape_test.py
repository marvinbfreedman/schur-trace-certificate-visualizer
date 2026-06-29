#!/usr/bin/env python3
import argparse

from analytic_mixed_derivative import phi_even_derivatives, source_st, w_values


def score_e(s: float, t: float, omega: float) -> float:
    fs, fps, _ = phi_even_derivatives(s)
    ft, fpt, _ = phi_even_derivatives(t)
    if fs == 0.0 or ft == 0.0:
        return 0.0
    h_s = -fps / fs
    h_t = -fpt / ft
    w0, w1, w2 = w_values(s + t, omega)
    return w2 - w1 * (h_s + h_t) + w0 * h_s * h_t


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=160)
    args = parser.parse_args()

    xs = [args.xmin + (args.xmax - args.xmin) * i / (args.n - 1) for i in range(args.n)]
    min_e = float("inf")
    min_d = float("inf")
    min_e_at = None
    min_d_at = None
    negative_count = 0
    for s in xs:
        for t in xs:
            e = score_e(s, t, args.omega)
            d = source_st(s, t, args.omega)
            if e < min_e:
                min_e = e
                min_e_at = (s, t)
            if d < min_d:
                min_d = d
                min_d_at = (s, t)
            if d < 0.0:
                negative_count += 1
    print(f"omega={args.omega:g} grid=[{args.xmin:g},{args.xmax:g}] n={args.n}")
    print(f"  min E={min_e:.12e} at s={min_e_at[0]:.6g}, t={min_e_at[1]:.6g}")
    print(f"  min D={min_d:.12e} at s={min_d_at[0]:.6g}, t={min_d_at[1]:.6g}")
    print(f"  negative D count={negative_count} / {args.n * args.n}")


if __name__ == "__main__":
    main()
