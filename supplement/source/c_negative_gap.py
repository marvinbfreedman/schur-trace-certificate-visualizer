#!/usr/bin/env python3
import argparse

from legendre_certificate import (
    complement,
    project_to_legendre,
    trial_coefficients,
    weighted_kernel,
)
from legendre_tail_decay import tail_norms

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("c_negative_gap.py requires numpy; run with python, not python3") from exc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--max-basis", type=int, default=64)
    parser.add_argument("--cuts", type=int, nargs="+", default=[16, 20, 24, 28, 32, 40])
    parser.add_argument("--quad", type=int, default=220)
    parser.add_argument("--laguerre", type=int, default=160)
    parser.add_argument("--trial", choices=("center_exp5", "s_exp1", "one_exp5"), default="center_exp5")
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    pts, wts, boundary_w, _ = weighted_kernel(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    cmat = project_to_legendre(boundary_w, pts, wts, args.max_basis, args.L)
    print(
        f"C negative gap kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"basis={args.max_basis} quad={args.quad} trial={args.trial}"
    )
    for cut in args.cuts:
        if cut >= args.max_basis:
            continue
        cblock = cmat[:cut, :cut]
        emat = trial_coefficients(pts, wts, cut, args.L, args.trial)
        rest = complement(emat, cut)
        vals = np.linalg.eigvalsh((rest.T @ cblock @ rest + rest.T @ cblock.T @ rest) / 2.0)
        neg = vals[vals < -args.tol]
        nonneg = vals[vals >= -args.tol]
        tail = tail_norms(cmat, cut)["outside_op"]
        first_nonneg = float(nonneg[0]) if len(nonneg) else float("nan")
        last_neg = float(neg[-1]) if len(neg) else float("nan")
        print(
            f"cut={cut}: neg={len(neg)} last_neg={last_neg:.12e} "
            f"first_nonneg={first_nonneg:.12e} tail_op={tail:.12e} "
            f"tail/first_gap={tail / max(first_nonneg, 1e-300):.3e}"
        )
        low = " ".join(f"{v:.4e}" for v in vals[: min(8, len(vals))])
        print(f"  low={low}")


if __name__ == "__main__":
    main()
