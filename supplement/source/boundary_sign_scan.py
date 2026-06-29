#!/usr/bin/env python3
import argparse
import math

from second_order_cosh_ibp import finite_modes, scalar_base_values

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "boundary_sign_scan.py requires numpy; run with PYTHONPATH=/private/tmp/codex_numpy"
    ) from exc


def signed_weight(s, t, omega, sign):
    u = 0.5 * (s + t)
    if sign == 0:
        return u
    return u * math.exp(sign * omega * u)


def boundary_value(x, y, omega, modes, sign):
    s = 2.0 * x
    t = 2.0 * y
    b_s, a_s, _, mu_s, _ = scalar_base_values(s, modes)
    b_t, a_t, _, mu_t, _ = scalar_base_values(t, modes)
    return b_s * b_t * signed_weight(s, t, omega, sign) * mu_s * mu_t / (a_s + a_t)


def matrix(xs, omega, modes, sign):
    out = np.zeros((len(xs), len(xs)), dtype=float)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs[: i + 1]):
            value = boundary_value(float(x), float(y), omega, modes, sign)
            out[i, j] = out[j, i] = value
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw1", "raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", type=float, default=0.49)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--tol", type=float, default=1e-10)
    args = parser.parse_args()

    modes = finite_modes(args.kind)
    xs = np.linspace(0.0, args.xmax, args.n)
    print(f"kind={args.kind} omega={args.omega:g} x=[0,{args.xmax:g}] n={args.n}")
    for sign in (1, -1, 0):
        mat = matrix(xs, args.omega, modes, sign)
        vals = np.linalg.eigvalsh((mat + mat.T) / 2.0)
        label = "plain" if sign == 0 else f"exp{sign:+d}"
        print(
            f"  {label}: neg={(vals < -args.tol).sum()} "
            f"zero={(abs(vals) <= args.tol).sum()} pos={(vals > args.tol).sum()} "
            f"min={vals[0]:.12e} max={vals[-1]:.12e}"
        )
    plus = matrix(xs, args.omega, modes, 1)
    minus = matrix(xs, args.omega, modes, -1)
    cosh = 0.5 * (plus + minus)
    vals = np.linalg.eigvalsh((cosh + cosh.T) / 2.0)
    print(
        f"  cosh half-sum: neg={(vals < -args.tol).sum()} "
        f"zero={(abs(vals) <= args.tol).sum()} pos={(vals > args.tol).sum()} "
        f"min={vals[0]:.12e} max={vals[-1]:.12e}"
    )


if __name__ == "__main__":
    main()
