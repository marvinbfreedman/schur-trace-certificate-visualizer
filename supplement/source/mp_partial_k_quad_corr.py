#!/usr/bin/env python3
import argparse
import mpmath as mp

from mp_partial_k_corr import parse_points, partial_phi, weights_for


def partial_k_quad(a, b, omega, weights):
    m = (a + b) / 2
    u = (a - b) / 2
    lower = abs(m)

    def integrand(r):
        y = lower + r
        return (
            mp.mpf("0.5")
            * y
            * mp.cosh(2 * omega * y)
            * partial_phi(y + u, weights)
            * partial_phi(y - u, weights)
        )

    # The large-|m| diagonal integrals live in a very thin boundary layer at r=0.
    return mp.quad(integrand, [0, mp.mpf("1e-5"), mp.mpf("5e-5"), mp.mpf("2e-4"),
                               mp.mpf("1e-3"), mp.mpf("5e-3"), mp.mpf("2e-2"),
                               mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("2.0"), mp.inf])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--points", required=True)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--vectors", action="store_true")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    xs = parse_points(args.points)
    weights = weights_for(args.kind)
    n = len(xs)
    mat = mp.matrix(n)
    diag = [partial_k_quad(x, x, omega, weights) for x in xs]
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            mat[i, j] = partial_k_quad(x, y, omega, weights) / mp.sqrt(diag[i] * diag[j])
    corr = (mat + mat.T) / 2
    if args.vectors:
        evals, evecs = mp.eigsy(corr, eigvals_only=False)
    else:
        evals = mp.eigsy(corr, eigvals_only=True)
        evecs = None
    print(f"quad kind={args.kind} omega={omega} n={n} dps={args.dps}")
    print("  alpha3=", mp.nstr(weights.get(3, mp.mpf(0)), 40))
    print("  diag min=", mp.nstr(min(diag), 30))
    print("  diag max=", mp.nstr(max(diag), 30))
    print("  corr min=", mp.nstr(evals[0], 40))
    print("  corr max=", mp.nstr(evals[-1], 30))
    print("  low=", " ".join(mp.nstr(value, 14) for value in evals[: min(8, n)]))
    if evecs is not None:
        print("  min eigenvector:")
        for x, coeff in zip(xs, [evecs[i, 0] for i in range(n)]):
            print("   ", mp.nstr(x, 14), mp.nstr(coeff, 18))


if __name__ == "__main__":
    main()
