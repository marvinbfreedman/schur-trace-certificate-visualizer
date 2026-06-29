#!/usr/bin/env python3
import argparse
import mpmath as mp

from mp_partial_k_corr import parse_points, weights_for


PIECES = ((4, mp.mpf("4.5")), (-6, mp.mpf("2.5")))


def mode_pieces(weights):
    pieces = []
    for n, weight in weights.items():
        n2 = n * n
        c = mp.pi * n2
        for base_coeff, lam in PIECES:
            if base_coeff == 4:
                coeff = weight * 4 * mp.pi * mp.pi * n2 * n2
            else:
                coeff = weight * -6 * mp.pi * n2
            pieces.append((coeff, lam, c))
    return pieces


def tail_moment(alpha, p):
    return alpha ** (-p) * mp.gammainc(p, alpha, mp.inf)


def exact_same_k(x, y, omega, pieces):
    # x,y >= 0, same-sign kernel K(x,y).
    total = mp.mpf(0)
    xy = x + y
    for coeff_i, lam_i, c_i in pieces:
        ex_i = mp.e ** (lam_i * x)
        cix = c_i * mp.e ** (2 * x)
        for coeff_j, lam_j, c_j in pieces:
            ex_j = mp.e ** (lam_j * y)
            cjy = c_j * mp.e ** (2 * y)
            alpha = cix + cjy
            for sigma in (-1, 1):
                p = (lam_i + lam_j + 2 * sigma * omega) / 2
                pref = (
                    coeff_i
                    * coeff_j
                    * ex_i
                    * ex_j
                    * mp.e ** (sigma * omega * xy)
                    / 16
                )
                i0 = tail_moment(alpha, p)
                ilog = mp.diff(lambda pp: tail_moment(alpha, pp), p)
                total += pref * (ilog + xy * i0)
    return total


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
    xs = [abs(x) for x in parse_points(args.points)]
    weights = weights_for(args.kind)
    pieces = mode_pieces(weights)
    n = len(xs)
    diag = [exact_same_k(x, x, omega, pieces) for x in xs]
    mat = mp.matrix(n)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            mat[i, j] = exact_same_k(x, y, omega, pieces) / mp.sqrt(diag[i] * diag[j])
    corr = (mat + mat.T) / 2
    if args.vectors:
        evals, evecs = mp.eigsy(corr, eigvals_only=False)
    else:
        evals = mp.eigsy(corr, eigvals_only=True)
        evecs = None
    print(f"exact same-sign kind={args.kind} omega={omega} n={n} dps={args.dps}")
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
