#!/usr/bin/env python3
import argparse

try:
    import mpmath as mp
except ImportError as exc:
    raise SystemExit("gaussian_mixture_test.py requires mpmath; run with PYTHONPATH=/private/tmp/codex_numpy") from exc


def phi_mp(t):
    x = abs(mp.mpf(t))
    total = mp.mpf("0")
    for n in range(1, 120):
        n2 = n * n
        term = (
            2
            * (2 * mp.pi * mp.pi * n2 * n2 * mp.e ** (mp.mpf("4.5") * x)
               - 3 * mp.pi * n2 * mp.e ** (mp.mpf("2.5") * x))
            * mp.e ** (-mp.pi * n2 * mp.e ** (2 * x))
        )
        total += term
        if n > 8 and abs(term) < mp.mpf("1e-70") * max(1, abs(total)):
            break
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", default="0.001,0.01,0.1,0.5,1,2")
    parser.add_argument("--max-order", type=int, default=8)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    mp.mp.dps = args.dps

    def f(x):
        return phi_mp(mp.sqrt(x))

    for text in args.points.split(","):
        x = mp.mpf(text)
        print(f"x={text}")
        for order in range(1, args.max_order + 1):
            signed = (-1) ** order * mp.diff(f, x, order)
            print(f"  order={order} signed={mp.nstr(signed, 20)}")


if __name__ == "__main__":
    main()
