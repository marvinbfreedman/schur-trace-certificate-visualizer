#!/usr/bin/env python3
import argparse
import mpmath as mp


def partial_phi(t, weights):
    x = abs(t)
    e2 = mp.e ** (2 * x)
    total = mp.mpf(0)
    for n, weight in weights.items():
        n2 = n * n
        c = mp.pi * n2
        total += weight * (
            4 * mp.pi * mp.pi * n2 * n2 * mp.e ** (mp.mpf("4.5") * x - c * e2)
            - 6 * mp.pi * n2 * mp.e ** (mp.mpf("2.5") * x - c * e2)
        )
    return total


def mode_right_derivative(n):
    n2 = n * n
    c = mp.pi * n2
    total = mp.mpf(0)
    for coeff, lam in (
        (4 * mp.pi * mp.pi * n2 * n2, mp.mpf("4.5")),
        (-6 * mp.pi * n2, mp.mpf("2.5")),
    ):
        term = coeff * mp.e ** (-c)
        total += term * (lam - 2 * c)
    return total


def weights_for(kind):
    if kind == "raw2":
        return {1: mp.mpf(1), 2: mp.mpf(1)}
    if kind == "raw3":
        return {1: mp.mpf(1), 2: mp.mpf(1), 3: mp.mpf(1)}
    if kind == "tilde3":
        alpha3 = -(mode_right_derivative(1) + mode_right_derivative(2)) / mode_right_derivative(3)
        return {1: mp.mpf(1), 2: mp.mpf(1), 3: alpha3}
    raise ValueError(kind)


def simpson_grid(rmax, intervals):
    if intervals % 2:
        raise ValueError("intervals must be even")
    h = rmax / intervals
    rs = [h * i for i in range(intervals + 1)]
    ws = []
    for i in range(intervals + 1):
        factor = 1 if i == 0 or i == intervals else 4 if i % 2 else 2
        ws.append(h * factor / 3)
    return rs, ws


def partial_k(a, b, omega, weights, rs, ws):
    m = (a + b) / 2
    u = (a - b) / 2
    lower = abs(m)
    total = mp.mpf(0)
    for r, w in zip(rs, ws):
        y = lower + r
        total += (
            mp.mpf("0.5")
            * w
            * y
            * mp.cosh(2 * omega * y)
            * partial_phi(y + u, weights)
            * partial_phi(y - u, weights)
        )
    return total


def parse_points(text):
    return [mp.mpf(piece) for piece in text.replace(",", " ").split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("raw2", "raw3", "tilde3"), default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--points", required=True)
    parser.add_argument("--rmax", default="9")
    parser.add_argument("--intervals", type=int, default=600)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--vectors", action="store_true")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    omega = mp.mpf(args.omega)
    xs = parse_points(args.points)
    weights = weights_for(args.kind)
    rs, ws = simpson_grid(mp.mpf(args.rmax), args.intervals)
    n = len(xs)
    mat = mp.matrix(n)
    diag = [partial_k(x, x, omega, weights, rs, ws) for x in xs]
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            denom = mp.sqrt(diag[i] * diag[j])
            mat[i, j] = partial_k(x, y, omega, weights, rs, ws) / denom
    corr = (mat + mat.T) / 2
    if args.vectors:
        evals, evecs = mp.eigsy(corr, eigvals_only=False)
    else:
        evals = mp.eigsy(corr, eigvals_only=True)
        evecs = None
    print(f"kind={args.kind} omega={omega} n={n} dps={args.dps}")
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
