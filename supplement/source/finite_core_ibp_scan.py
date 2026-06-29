#!/usr/bin/env python3
import argparse
import math

from mp_partial_k_corr import weights_for

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("finite_core_ibp_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


PI = math.pi


def pieces(kind: str, omega: float, eps: int):
    # F_eps(v) = exp(eps*omega*v/2) phi(v/2)
    #          = sum A_i exp(beta_i v - c_i exp(v)).
    weights = weights_for(kind)
    out = []
    for n, weight in weights.items():
        wf = float(weight)
        n2 = n * n
        c = PI * n2
        for coeff, lam in (
            (4.0 * PI * PI * n2 * n2, 4.5),
            (-6.0 * PI * n2, 2.5),
        ):
            beta = 0.5 * (lam + eps * omega)
            out.append((wf * coeff, beta, c))
    return out


def f_derivs(v: float, pcs):
    ev = math.exp(v)
    f0 = 0.0
    f1 = 0.0
    f2 = 0.0
    for coeff, beta, c in pcs:
        term = coeff * math.exp(beta * v - c * ev)
        ld = beta - c * ev
        f0 += term
        f1 += term * ld
        f2 += term * (ld * ld - c * ev)
    return f0, f1, f2


def score_and_deriv(v: float, pcs):
    f0, f1, f2 = f_derivs(v, pcs)
    a = -f1 / f0
    ap = -(f2 * f0 - f1 * f1) / (f0 * f0)
    return a, ap


def ibp_coeff(s: float, t: float, pcs):
    # C=(s+t)/(2(a(s)+a(t))); D=(partial_s+partial_t)C.
    a, ap = score_and_deriv(s, pcs)
    b, bp = score_and_deriv(t, pcs)
    q = a + b
    return (q - 0.5 * (s + t) * (ap + bp)) / (q * q)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--eps", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=8.0)
    parser.add_argument("--n", type=int, default=80)
    args = parser.parse_args()

    pcs = pieces(args.kind, args.omega, args.eps)
    xs = np.linspace(args.vmin, args.vmax, args.n)
    min_coeff = float("inf")
    min_at = None
    mat = np.zeros((args.n, args.n), dtype=float)
    for i, s in enumerate(xs):
        for j, t in enumerate(xs[: i + 1]):
            value = ibp_coeff(float(s), float(t), pcs)
            mat[i, j] = mat[j, i] = value
            if value < min_coeff:
                min_coeff = value
                min_at = (float(s), float(t))
    evals = np.linalg.eigvalsh(mat)
    print(
        f"kind={args.kind} omega={args.omega:g} eps={args.eps:+d} "
        f"v=[{args.vmin:g},{args.vmax:g}] n={args.n}"
    )
    print(f"  min coeff={min_coeff:.12e} at s={min_at[0]:.6g} t={min_at[1]:.6g}")
    print(f"  coeff-kernel min eig={evals[0]:.12e} max eig={evals[-1]:.12e}")
    for v in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0):
        if args.vmin <= v <= args.vmax:
            a, ap = score_and_deriv(v, pcs)
            print(f"  score v={v:g}: a={a:.12e} a'={ap:.12e} rho={v / a if a else float('nan'):.12e}")


if __name__ == "__main__":
    main()
