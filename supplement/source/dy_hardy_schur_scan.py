#!/usr/bin/env python3
import argparse
import math

from dy_hardy_split import build_matrices, eigvals, generalized_ratios, schur_relative_to_negative


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-min", type=float, default=1.25)
    parser.add_argument("--p-max", type=float, default=2.5)
    parser.add_argument("--p-count", type=int, default=26)
    parser.add_argument("--c", type=float, default=math.pi)
    parser.add_argument("--L", type=float, default=12.0)
    parser.add_argument("--s-order", type=int, default=80)
    parser.add_argument("--umax", type=float, default=16.0)
    parser.add_argument("--u-order", type=int, default=320)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    print(
        f"Dy Hardy Schur scan p=[{args.p_min:g},{args.p_max:g}] "
        f"count={args.p_count} L={args.L:g} s_order={args.s_order} "
        f"umax={args.umax:g} u_order={args.u_order}"
    )
    worst_ratio = (float("inf"), None)
    worst_schur = (float("inf"), None)
    worst_full = (float("inf"), None)
    for i in range(args.p_count):
        if args.p_count == 1:
            p = args.p_min
        else:
            p = args.p_min + (args.p_max - args.p_min) * i / (args.p_count - 1)
        mats = build_matrices(p, args.c, args.L, args.s_order, args.umax, args.u_order)
        lower_vals = eigvals(mats["lower"])
        hardy_vals = eigvals(mats["P_plus_R"])
        full_vals = eigvals(mats["full_expanded"])
        _, ratios = generalized_ratios(4.0 * mats["P_plus_R"], mats["lower"], args.tol)
        ratio_min = float("nan") if ratios is None else float(ratios[0])
        schur, defect, rank = schur_relative_to_negative(
            mats["full_expanded"], mats["lower"], args.tol
        )
        schur_min = float("nan")
        if schur is not None:
            schur_min = float(eigvals(schur)[0])
        if ratios is not None and ratio_min < worst_ratio[0]:
            worst_ratio = (ratio_min, p)
        if not math.isnan(schur_min) and schur_min < worst_schur[0]:
            worst_schur = (schur_min, p)
        if full_vals[0] < worst_full[0]:
            worst_full = (float(full_vals[0]), p)
        print(
            f"  p={p:.6f} "
            f"Hmin={hardy_vals[0]: .3e} "
            f"Lmin={lower_vals[0]: .3e} Lneg={(lower_vals < -args.tol).sum():2d} "
            f"Emin={full_vals[0]: .3e} "
            f"ratio={ratio_min: .3e} "
            f"Schur={schur_min: .3e} defect={defect: .1e} rank={rank}"
        )
    print(
        f"worst ratio={worst_ratio[0]:.12e} at p={worst_ratio[1]:.6g}; "
        f"worst Schur={worst_schur[0]:.12e} at p={worst_schur[1]:.6g}; "
        f"worst Emin={worst_full[0]:.12e} at p={worst_full[1]:.6g}"
    )


if __name__ == "__main__":
    main()
