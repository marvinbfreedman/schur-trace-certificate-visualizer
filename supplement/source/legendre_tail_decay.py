#!/usr/bin/env python3
import argparse

from legendre_certificate import project_to_legendre, weighted_kernel

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "legendre_tail_decay.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def tail_norms(mat, cut):
    tail = mat.copy()
    tail[:cut, :cut] = 0.0
    high = mat[cut:, cut:]
    cross = np.zeros_like(mat)
    cross[:cut, cut:] = mat[:cut, cut:]
    cross[cut:, :cut] = mat[cut:, :cut]
    return {
        "outside_frob": float(np.linalg.norm(tail, "fro")),
        "outside_op": float(np.linalg.norm(tail, 2)),
        "high_frob": float(np.linalg.norm(high, "fro")),
        "high_op": float(np.linalg.norm(high, 2)) if high.size else 0.0,
        "cross_frob": float(np.linalg.norm(cross, "fro")),
        "cross_op": float(np.linalg.norm(cross, 2)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--L", type=float, default=8.0)
    parser.add_argument("--max-basis", type=int, default=60)
    parser.add_argument("--cuts", type=int, nargs="+", default=[12, 16, 20, 24, 28, 32, 40])
    parser.add_argument("--quad", type=int, default=260)
    parser.add_argument("--laguerre", type=int, default=180)
    args = parser.parse_args()

    pts, wts, boundary_w, full_w = weighted_kernel(
        args.kind, args.omega, args.L, args.quad, args.laguerre
    )
    cmat = project_to_legendre(boundary_w, pts, wts, args.max_basis, args.L)
    kmat = project_to_legendre(full_w, pts, wts, args.max_basis, args.L)
    tmat = kmat - cmat
    print(
        f"legendre tail kind={args.kind} omega={args.omega:g} L={args.L:g} "
        f"basis={args.max_basis} quad={args.quad} laguerre={args.laguerre}"
    )
    for name, mat in (("K", kmat), ("C", cmat), ("T", tmat)):
        print(f"{name}: total_frob={np.linalg.norm(mat, 'fro'):.12e} op={np.linalg.norm(mat, 2):.12e}")
        for cut in args.cuts:
            if cut >= args.max_basis:
                continue
            data = tail_norms(mat, cut)
            print(
                f"  cut={cut}: outside_op={data['outside_op']:.12e} "
                f"outside_frob={data['outside_frob']:.12e} "
                f"cross_op={data['cross_op']:.12e} high_op={data['high_op']:.12e}"
            )


if __name__ == "__main__":
    main()
