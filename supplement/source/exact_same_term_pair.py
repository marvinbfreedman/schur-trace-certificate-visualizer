#!/usr/bin/env python3
import argparse
import mpmath as mp

from exact_same_k_corr import exact_same_k, mode_pieces
from mp_partial_k_corr import parse_points


def pair_weights(n, m, sym):
    if sym and n != m:
        # The kernel formula is bilinear in phi. Use a marker handled below.
        return None
    return {n: mp.mpf(1)} if n == m else {n: mp.mpf(1), m: mp.mpf(0)}


def exact_pair_k(x, y, omega, n, m, sym):
    pieces_n = mode_pieces({n: mp.mpf(1)})
    pieces_m = mode_pieces({m: mp.mpf(1)})
    # exact_same_k expects one list for both sides; duplicate the loop here.
    def one_side(left_pieces, right_pieces):
        total = mp.mpf(0)
        xy = x + y
        for coeff_i, lam_i, c_i in left_pieces:
            ex_i = mp.e ** (lam_i * x)
            cix = c_i * mp.e ** (2 * x)
            for coeff_j, lam_j, c_j in right_pieces:
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
                    i0 = alpha ** (-p) * mp.gammainc(p, alpha, mp.inf)
                    ilog = mp.diff(lambda pp: alpha ** (-pp) * mp.gammainc(pp, alpha, mp.inf), p)
                    total += pref * (ilog + xy * i0)
        return total

    value = one_side(pieces_n, pieces_m)
    if sym and n != m:
        value += one_side(pieces_m, pieces_n)
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--pairs", default="1,1;1,2;2,2;1,3;2,3;3,3")
    parser.add_argument("--points", default="0 0.1 0.2 0.5 1 1.5 2 2.5")
    parser.add_argument("--dps", type=int, default=70)
    args = parser.parse_args()

    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    xs = [abs(x) for x in parse_points(args.points)]
    for spec in args.pairs.split(";"):
        n_text, m_text = spec.split(",")
        n = int(n_text)
        m = int(m_text)
        size = len(xs)
        mat = mp.matrix(size)
        diag = []
        # Use raw pair diagonal normalization, not meaningful if diagonal changes sign.
        for x in xs:
            diag.append(exact_pair_k(x, x, omega, n, m, True))
        for i, x in enumerate(xs):
            for j, y in enumerate(xs):
                mat[i, j] = exact_pair_k(x, y, omega, n, m, True)
        evals = mp.eigsy((mat + mat.T) / 2, eigvals_only=True)
        print(f"pair ({n},{m}) omega={omega}: min={mp.nstr(evals[0], 18)} max={mp.nstr(evals[-1], 18)}")
        print(f"  diag min={mp.nstr(min(diag), 18)} max={mp.nstr(max(diag), 18)}")


if __name__ == "__main__":
    main()
