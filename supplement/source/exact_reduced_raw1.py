#!/usr/bin/env python3
import argparse
import mpmath as mp


def points_for(args):
    if args.grid == "linear":
        return [args.vmin + (args.vmax - args.vmin) * i / (args.n - 1) for i in range(args.n)]
    if args.grid == "quadratic":
        return [
            args.vmin + (args.vmax - args.vmin) * (mp.mpf(i) / (args.n - 1)) ** 2
            for i in range(args.n)
        ]
    if args.grid == "geometric":
        if args.vmin <= 0:
            raise ValueError("geometric grid requires --vmin > 0")
        ratio = args.vmax / args.vmin
        return [args.vmin * ratio ** (mp.mpf(i) / (args.n - 1)) for i in range(args.n)]
    raise ValueError(args.grid)


def mode_mu(v):
    c = mp.pi
    ev = mp.e**v
    return 4 * c * c * ev * ev - 6 * c * ev


def mode_score(v):
    return mp.pi * mp.e**v - mp.mpf("0.25")


def boundary_red(s, t, omega):
    u = (s + t) / 2
    if u == 0:
        return mp.mpf(0)
    return u * mp.cosh(omega * u) / (mode_score(s) + mode_score(t))


def tail_moment(alpha, p):
    # int_1^inf q^p exp(-alpha q) dq
    return alpha ** (-(p + 1)) * mp.gammainc(p + 1, alpha, mp.inf)


def exact_full_red(s, t, omega):
    c = mp.pi
    es = mp.e**s
    et = mp.e**t
    alpha = c * (es + et)
    center = (s + t) / 2
    mu_s = mode_mu(s)
    mu_t = mode_mu(t)
    # mu(s+u)/mu(s) = b1(s) q + b2(s) q^2.
    b_s = {
        1: -6 * c * es / mu_s,
        2: 4 * c * c * es * es / mu_s,
    }
    b_t = {
        1: -6 * c * et / mu_t,
        2: 4 * c * c * et * et / mu_t,
    }
    total = mp.mpf(0)
    for i, coeff_i in b_s.items():
        for j, coeff_j in b_t.items():
            coeff = coeff_i * coeff_j
            for sigma in (-1, 1):
                p = i + j - mp.mpf("0.5") + sigma * omega
                i0 = tail_moment(alpha, p)
                ilog = mp.diff(lambda pp: tail_moment(alpha, pp), p)
                total += (
                    mp.mpf("0.5")
                    * coeff
                    * mp.e ** (alpha + sigma * omega * center)
                    * (center * i0 + ilog)
                )
    return total


def eigvals(mat):
    return mp.eigsy((mat + mat.T) / 2, eigvals_only=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--vmin", default="0")
    parser.add_argument("--vmax", default="5.2")
    parser.add_argument("--n", type=int, default=30)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--tol", default="1e-40")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    args.vmin = mp.mpf(args.vmin)
    args.vmax = mp.mpf(args.vmax)
    tol = mp.mpf(args.tol)
    pts = points_for(args)

    boundary = mp.matrix(args.n)
    full = mp.matrix(args.n)
    for r, s in enumerate(pts):
        for c, t in enumerate(pts):
            boundary[r, c] = boundary_red(s, t, omega)
            full[r, c] = exact_full_red(s, t, omega)
    tail = full - boundary
    print(
        f"exact reduced raw1 omega={omega} grid={args.grid} "
        f"v=[{args.vmin},{args.vmax}] n={args.n} dps={args.dps}"
    )
    for name, mat in (("boundary", boundary), ("tail", tail), ("full", full)):
        vals = eigvals(mat)
        vals_list = [vals[i] for i in range(args.n)]
        neg = sum(1 for value in vals_list if value < -tol)
        zero = sum(1 for value in vals_list if abs(value) <= tol)
        pos = sum(1 for value in vals_list if value > tol)
        print(f"  {name}: neg={neg} zero={zero} pos={pos}")
        print("    low=", " ".join(mp.nstr(value, 18) for value in vals_list[: min(10, args.n)]))
        print("    high=", " ".join(mp.nstr(value, 18) for value in vals_list[-min(6, args.n):]))


if __name__ == "__main__":
    main()
