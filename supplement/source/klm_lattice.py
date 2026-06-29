#!/usr/bin/env python3
import argparse
import math
import random
import time

from klm_test import SymplecticFourier, hermitian_real_form, jacobi_eigenvalues_symmetric

try:
    import numpy as np
except ImportError:
    np = None


def lattice_points(rows: int, cols: int, dy: float, deta: float):
    cy = 0.5 * (rows - 1)
    ce = 0.5 * (cols - 1)
    return [((r - cy) * dy, (c - ce) * deta) for r in range(rows) for c in range(cols)]


def precompute_kernel(points, qfun, hbar: float):
    n = len(points)
    kernel = [[0j for _ in range(n)] for _ in range(n)]
    for j, (yj, etaj) in enumerate(points):
        for k, (yk, etak) in enumerate(points):
            q = qfun(yj - yk, etaj - etak)
            symp = etaj * yk - yj * etak
            kernel[j][k] = q * complex(math.cos(0.5 * hbar * symp), math.sin(0.5 * hbar * symp))
    return kernel


def matvec(kernel, vector):
    out = []
    for row in kernel:
        total = 0j
        for a, x in zip(row, vector):
            total += a * x
        out.append(total)
    return out


def dotc(a, b):
    return sum(x.conjugate() * y for x, y in zip(a, b))


def norm(v):
    return math.sqrt(max(0.0, dotc(v, v).real))


def normalize(v):
    nrm = norm(v)
    if nrm == 0.0:
        raise ValueError("zero vector")
    return [x / nrm for x in v]


def lanczos_min(kernel, steps: int, seed: int):
    n = len(kernel)
    rng = random.Random(seed)
    q = normalize([complex(rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)) for _ in range(n)])
    basis = []
    alphas = []
    betas = []
    prev = [0j for _ in range(n)]
    beta_prev = 0.0

    for _ in range(min(steps, n)):
        basis.append(q)
        z = matvec(kernel, q)
        alpha = dotc(q, z).real
        z = [zi - alpha * qi - beta_prev * pi for zi, qi, pi in zip(z, q, prev)]

        # Full reorthogonalization keeps the small negative-witness search honest.
        for b in basis:
            coeff = dotc(b, z)
            z = [zi - coeff * bi for zi, bi in zip(z, b)]

        beta = norm(z)
        alphas.append(alpha)
        if beta < 1e-13:
            break
        betas.append(beta)
        prev = q
        q = [zi / beta for zi in z]
        beta_prev = beta

    m = len(alphas)
    tri = [[0.0 for _ in range(m)] for _ in range(m)]
    for i, alpha in enumerate(alphas):
        tri[i][i] = alpha
    for i, beta in enumerate(betas[: m - 1]):
        tri[i][i + 1] = beta
        tri[i + 1][i] = beta
    return jacobi_eigenvalues_symmetric(tri)[0], m


def exact_min(kernel):
    if np is not None:
        return float(np.linalg.eigvalsh(np.array(kernel, dtype=np.complex128))[0]), len(kernel)
    eigs = jacobi_eigenvalues_symmetric(hermitian_real_form(kernel), tol=1e-12)
    return eigs[0], len(kernel)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--hbar", type=float, default=1.0)
    parser.add_argument("--rows", type=int, default=7)
    parser.add_argument("--cols", type=int, default=7)
    parser.add_argument("--dy-values", default="0.25,0.5,0.75,1,1.25,1.5")
    parser.add_argument("--deta-values", default="0.25,0.5,0.75,1,1.25,1.5,2,2.5,3")
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--ymax", type=float, default=6.5)
    parser.add_argument("--lanczos-steps", type=int, default=60)
    parser.add_argument("--seed", type=int, default=20260528)
    parser.add_argument("--solver", choices=("exact", "lanczos"), default="exact")
    args = parser.parse_args()

    start = time.time()
    qfun = SymplecticFourier(args.omega, args.ymax, args.intervals)
    dy_values = [float(x) for x in args.dy_values.split(",") if x]
    deta_values = [float(x) for x in args.deta_values.split(",") if x]
    best = (float("inf"), None, None, None)

    backend = "numpy" if np is not None else "pure-python"
    print(f"omega={args.omega:g} hbar={args.hbar:g} q00={qfun.q00:.12e} backend={backend}")
    for dy in dy_values:
        for deta in deta_values:
            points = lattice_points(args.rows, args.cols, dy, deta)
            kernel = precompute_kernel(points, qfun, args.hbar)
            if args.solver == "exact":
                lam, used = exact_min(kernel)
            else:
                lam, used = lanczos_min(kernel, args.lanczos_steps, args.seed)
            print(
                f"  rows={args.rows} cols={args.cols} dy={dy:g} deta={deta:g} "
                f"lambda_min~{lam:.12e} solver={args.solver}:{used}"
            )
            if lam < best[0]:
                best = (lam, dy, deta, used)

    print(
        f"best lambda_min~{best[0]:.12e} dy={best[1]:g} deta={best[2]:g} "
        f"solver={args.solver}:{best[3]} q-cache={len(qfun.cache)} elapsed={time.time() - start:.2f}s"
    )


if __name__ == "__main__":
    main()
