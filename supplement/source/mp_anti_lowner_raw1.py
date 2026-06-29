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


def raw1_kernel(s, t):
    if s == 0 and t == 0:
        return mp.mpf(0)
    return (s + t) / (mp.pi * (mp.e**s + mp.e**t) - mp.mpf("0.5"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vmin", default="0")
    parser.add_argument("--vmax", default="5.2")
    parser.add_argument("--n", type=int, default=30)
    parser.add_argument("--grid", choices=("linear", "quadratic", "geometric"), default="linear")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--tol", default="1e-40")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    args.vmin = mp.mpf(args.vmin)
    args.vmax = mp.mpf(args.vmax)
    tol = mp.mpf(args.tol)

    pts = points_for(args)
    mat = mp.matrix(args.n)
    for i, s in enumerate(pts):
        for j, t in enumerate(pts):
            mat[i, j] = raw1_kernel(s, t)
    evals = mp.eigsy((mat + mat.T) / 2, eigvals_only=True)
    neg = sum(1 for value in evals if value < -tol)
    zero = sum(1 for value in evals if abs(value) <= tol)
    pos = sum(1 for value in evals if value > tol)
    print(
        f"raw1 mp grid={args.grid} v=[{args.vmin},{args.vmax}] "
        f"n={args.n} dps={args.dps} tol={tol}"
    )
    print(f"  inertia: neg={neg} zero={zero} pos={pos}")
    print("  low=", " ".join(mp.nstr(value, 18) for value in evals[: min(12, args.n)]))
    eval_list = [evals[i] for i in range(args.n)]
    print("  high=", " ".join(mp.nstr(value, 18) for value in eval_list[-min(8, args.n):]))


if __name__ == "__main__":
    main()
