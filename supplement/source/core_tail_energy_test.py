#!/usr/bin/env python3
import argparse

from partial_phi_mixed_test import full_partial_h
from ranged_phi_mixed_test import full_ranged_h
from klm_test import simpson_grid

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("core_tail_energy_test.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def fill_matrix(xs, omega: float, rs, ws, kind: str, nmax: int):
    mat = np.zeros((len(xs), len(xs)), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            if kind == "partial":
                value = full_partial_h(float(x), float(y), omega, nmax, rs, ws)
            elif kind == "tail":
                value = full_ranged_h(float(x), float(y), omega, nmax + 1, 12, rs, ws)
            elif kind == "full":
                value = full_partial_h(float(x), float(y), omega, 12, rs, ws)
            else:
                raise ValueError(kind)
            mat[i, j] = mat[j, i] = value
    return mat


def generalized_norm(core, perturbation, relative_cutoff: float):
    evals, basis = np.linalg.eigh(core)
    keep = evals > relative_cutoff * max(evals[-1], 1.0)
    inv_sqrt = basis[:, keep] / np.sqrt(evals[keep])
    reduced = inv_sqrt.T @ perturbation @ inv_sqrt
    spectrum = np.linalg.eigvalsh(reduced)
    null_basis = basis[:, ~keep]
    if null_basis.shape[1]:
        null_spectrum = np.linalg.eigvalsh(null_basis.T @ perturbation @ null_basis)
    else:
        null_spectrum = np.array([])
    cross_norm = 0.0
    if null_basis.shape[1] and keep.sum():
        cross_norm = float(np.linalg.svd(null_basis.T @ perturbation @ inv_sqrt, compute_uv=False)[0])
    return evals, keep, spectrum, null_spectrum, cross_norm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--n-grid", type=int, default=64)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--core-nmax", type=int, default=2)
    parser.add_argument("--rmax", type=float, default=8.0)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--relative-cutoff", type=float, default=1e-12)
    args = parser.parse_args()

    xs = np.linspace(-args.xmax, args.xmax, args.n_grid)
    rs, ws = simpson_grid(args.rmax, args.intervals)
    core = fill_matrix(xs, args.omega, rs, ws, "partial", args.core_nmax)
    tail = fill_matrix(xs, args.omega, rs, ws, "tail", args.core_nmax)
    full = core + tail

    core_evals, keep, tail_spectrum, null_spectrum, cross_norm = generalized_norm(
        core, tail, args.relative_cutoff
    )
    full_evals = np.linalg.eigvalsh(full)
    tail_evals = np.linalg.eigvalsh(tail)

    print(f"omega={args.omega:g} grid={args.n_grid} x=[{-args.xmax:g},{args.xmax:g}] core<= {args.core_nmax}")
    print(f"  core min={core_evals[0]:.12e} max={core_evals[-1]:.12e}")
    print(f"  tail min={tail_evals[0]:.12e} max={tail_evals[-1]:.12e}")
    print(f"  full min={full_evals[0]:.12e} max={full_evals[-1]:.12e}")
    print(f"  retained core rank={int(keep.sum())} / {args.n_grid}")
    print(f"  tail relative min={tail_spectrum[0]:.12e} max={tail_spectrum[-1]:.12e}")
    if len(null_spectrum):
        print(f"  tail on core-null min={null_spectrum[0]:.12e} max={null_spectrum[-1]:.12e}")
        print(f"  tail cross null/range norm={cross_norm:.12e}")


if __name__ == "__main__":
    main()
