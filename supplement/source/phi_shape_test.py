#!/usr/bin/env python3
import argparse

from klm_test import phi


def d1(x: float, h: float) -> float:
    return (phi(x + h) - phi(x - h)) / (2.0 * h)


def d2_log(x: float, h: float) -> float:
    import math

    return (math.log(phi(x + h)) - 2.0 * math.log(phi(x)) + math.log(phi(x - h))) / (h * h)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xmin", type=float, default=1e-4)
    parser.add_argument("--xmax", type=float, default=2.6)
    parser.add_argument("--n", type=int, default=400)
    parser.add_argument("--h", type=float, default=1e-5)
    args = parser.parse_args()

    xs = [args.xmin + (args.xmax - args.xmin) * i / (args.n - 1) for i in range(args.n)]
    scores = []
    ratios = []
    curvatures = []
    for x in xs:
        p = phi(x)
        score = -d1(x, args.h) / p
        scores.append(score)
        ratios.append(x / score)
        curvatures.append(-d2_log(x, args.h))

    print(f"x range [{xs[0]:.6g}, {xs[-1]:.6g}], n={len(xs)}")
    print(f"min Phi={min(phi(x) for x in xs):.12e}")
    print(f"min score-2x={min(score - 2.0 * x for score, x in zip(scores, xs)):.12e}")
    print(f"ratio x/score max={max(ratios):.12e}, min={min(ratios):.12e}")
    print(f"max diff ratio={max(ratios[i + 1] - ratios[i] for i in range(len(ratios) - 1)):.12e}")
    print(f"min -log(Phi)''={min(curvatures):.12e}")


if __name__ == "__main__":
    main()
