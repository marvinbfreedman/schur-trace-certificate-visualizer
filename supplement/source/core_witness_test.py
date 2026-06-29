#!/usr/bin/env python3
import argparse

from partial_phi_mixed_test import full_partial_h
from ranged_phi_mixed_test import full_ranged_h
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("core_witness_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def matrix_for(xs, omega: float, rs, ws, kind: str):
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            if kind == "core2":
                value = full_partial_h(float(x), float(y), omega, 2, rs, ws)
            elif kind == "core3":
                value = full_partial_h(float(x), float(y), omega, 3, rs, ws)
            elif kind == "full12":
                value = full_partial_h(float(x), float(y), omega, 12, rs, ws)
            else:
                raise ValueError(kind)
            mat[i, j] = mat[j, i] = value
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n-grid", type=int, default=64)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--rmax", type=float, default=9.0)
    parser.add_argument("--intervals", type=int, default=1000)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    core2 = matrix_for(xs, args.omega, rs, ws, "core2")
    evals, basis = np.linalg.eigh(core2)
    witness = basis[:, 0]

    core3 = matrix_for(xs, args.omega, rs, ws, "core3")
    full12 = matrix_for(xs, args.omega, rs, ws, "full12")
    correction3 = core3 - core2
    tail4 = full12 - core3

    def qform(mat):
        return float(witness @ mat @ witness)

    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}]")
    print(f"  lambda_min core2={evals[0]:.12e}")
    print(f"  witness core2={qform(core2):.12e}")
    print(f"  witness true n=3 correction={qform(correction3):.12e}")
    print(f"  witness core3={qform(core3):.12e}")
    print(f"  witness tail4={qform(tail4):.12e}")
    print("  largest witness entries:")
    order = np.argsort(-np.abs(witness))[:8]
    for idx in order:
        print(f"    x={xs[idx]: .12e} coeff={witness[idx]: .12e}")


if __name__ == "__main__":
    main()
