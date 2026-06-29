#!/usr/bin/env python3
import argparse
import cmath
import math
import random
import time

PI = math.pi


def phi(t: float) -> float:
    x = abs(t)
    if x >= 8.0:
        return 0.0

    ex2 = math.exp(2.0 * x)
    total = 0.0
    for n in range(1, 80):
        n2 = n * n
        term = (
            2.0
            * (2.0 * PI * PI * n2 * n2 * math.exp(4.5 * x) - 3.0 * PI * n2 * math.exp(2.5 * x))
            * math.exp(-PI * n2 * ex2)
        )
        total += term
        if n > 4 and abs(term) < 1e-16 * max(1.0, abs(total)):
            break
    return total


def simpson_grid(ymax: float, intervals: int):
    if intervals % 2:
        intervals += 1
    h = ymax / intervals
    ys = [i * h for i in range(intervals + 1)]
    ws = []
    for i in range(intervals + 1):
        if i == 0 or i == intervals:
            w = 1.0
        elif i % 2:
            w = 4.0
        else:
            w = 2.0
        ws.append(w * h / 3.0)
    return ys, ws


class SymplecticFourier:
    """Q(y, eta) = F_sigma[sigma_omega](y, eta), normalized by Q(0, 0)."""

    def __init__(self, omega: float, ymax: float, intervals: int):
        self.omega = omega
        self.ys, self.ws = simpson_grid(ymax, intervals)
        self.cache = {}
        self.q00 = self.raw(0.0, 0.0)

    def raw(self, ydisp: float, eta: float) -> float:
        a = abs(ydisp)
        b = abs(eta)
        key = (round(a, 11), round(b, 11))
        if key in self.cache:
            return self.cache[key]

        half = 0.5 * a
        total = 0.0
        for y, w in zip(self.ys, self.ws):
            if b < 1e-10:
                osc = y
            else:
                osc = math.sin(b * y) / b
            total += w * y * math.cosh(2.0 * self.omega * y) * phi(y + half) * phi(y - half) * osc

        value = 2.0 * PI * total
        self.cache[key] = value
        return value

    def __call__(self, ydisp: float, eta: float) -> float:
        return self.raw(ydisp, eta) / self.q00


def klm_matrix(points, qfun, hbar: float):
    n = len(points)
    mat = [[0j for _ in range(n)] for _ in range(n)]
    for j, (yj, etaj) in enumerate(points):
        for k, (yk, etak) in enumerate(points):
            q = qfun(yj - yk, etaj - etak)
            symp = etaj * yk - yj * etak
            mat[j][k] = q * cmath.exp(0.5j * hbar * symp)
    return mat


def hermitian_real_form(mat):
    n = len(mat)
    out = [[0.0 for _ in range(2 * n)] for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            z = mat[i][j]
            out[i][j] = z.real
            out[i][j + n] = -z.imag
            out[i + n][j] = z.imag
            out[i + n][j + n] = z.real
    return out


def jacobi_eigenvalues_symmetric(a, tol=1e-11, max_sweeps=None):
    n = len(a)
    if max_sweeps is None:
        max_sweeps = max(80, 10 * n * n)
    for _ in range(max_sweeps):
        p, q = 0, 1
        max_off = 0.0
        for i in range(n):
            row = a[i]
            for j in range(i + 1, n):
                v = abs(row[j])
                if v > max_off:
                    max_off = v
                    p, q = i, j
        if max_off < tol:
            break

        app = a[p][p]
        aqq = a[q][q]
        apq = a[p][q]
        if abs(apq) < tol:
            continue
        tau = (aqq - app) / (2.0 * apq)
        sign = 1.0 if tau >= 0.0 else -1.0
        t = sign / (abs(tau) + math.sqrt(1.0 + tau * tau))
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c

        for k in range(n):
            if k == p or k == q:
                continue
            akp = a[k][p]
            akq = a[k][q]
            new_kp = c * akp - s * akq
            new_kq = s * akp + c * akq
            a[k][p] = new_kp
            a[p][k] = new_kp
            a[k][q] = new_kq
            a[q][k] = new_kq

        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        a[p][q] = 0.0
        a[q][p] = 0.0

    return sorted(a[i][i] for i in range(n))


def min_klm_eigenvalue(points, qfun, hbar: float):
    eigs = jacobi_eigenvalues_symmetric(hermitian_real_form(klm_matrix(points, qfun, hbar)))
    return eigs[0]


def rectangle_points(a: float, b: float):
    return [(0.0, 0.0), (a, 0.0), (0.0, b), (a, b)]


def triangle_points(a: float, b: float):
    return [(0.0, 0.0), (a, 0.0), (0.0, b)]


def structured_search(qfun, hbar: float, step: float, limit: float):
    best = (float("inf"), None, None)
    count = int(round(limit / step))
    for ia in range(1, count + 1):
        a = ia * step
        for ib in range(1, count + 1):
            b = ib * step
            for name, points in (("tri", triangle_points(a, b)), ("rect", rectangle_points(a, b))):
                lam = min_klm_eigenvalue(points, qfun, hbar)
                if lam < best[0]:
                    best = (lam, name, points)
    return best


def separated(points, min_sep: float) -> bool:
    if min_sep <= 0.0:
        return True
    min_sep2 = min_sep * min_sep
    for i, (yi, etai) in enumerate(points):
        for yj, etaj in points[:i]:
            dy = yi - yj
            de = etai - etaj
            if dy * dy + de * de < min_sep2:
                return False
    return True


def random_search(qfun, hbar: float, samples: int, npoints: int, radius: float, seed: int, min_sep: float):
    rng = random.Random(seed)
    best = (float("inf"), None)
    tried = 0
    accepted = 0
    while accepted < samples and tried < max(1000, 100 * samples):
        tried += 1
        points = [(rng.uniform(-radius, radius), rng.uniform(-radius, radius)) for _ in range(npoints)]
        if not separated(points, min_sep):
            continue
        accepted += 1
        lam = min_klm_eigenvalue(points, qfun, hbar)
        if lam < best[0]:
            best = (lam, points)
    return best, accepted


def fmt_points(points):
    return "[" + ", ".join(f"({y:.6g}, {eta:.6g})" for y, eta in points) + "]"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omegas", default="0,0.25,0.49")
    parser.add_argument("--ymax", type=float, default=6.5)
    parser.add_argument("--intervals", type=int, default=700)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--limit", type=float, default=4.0)
    parser.add_argument("--random-samples", type=int, default=80)
    parser.add_argument("--random-points", type=int, default=6)
    parser.add_argument("--random-radius", type=float, default=4.0)
    parser.add_argument("--min-sep", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20260520)
    parser.add_argument("--hbar", type=float, default=1.0)
    args = parser.parse_args()

    for omega_text in args.omegas.split(","):
        omega = float(omega_text)
        start = time.time()
        qfun = SymplecticFourier(omega, args.ymax, args.intervals)
        print(f"omega={omega:g} q00={qfun.q00:.12e}")

        lam, name, points = structured_search(qfun, args.hbar, args.step, args.limit)
        print(f"  structured best: lambda_min={lam:.12e} type={name} points={fmt_points(points)}")

        if args.random_samples > 0:
            (lam, points), accepted = random_search(
                qfun,
                args.hbar,
                args.random_samples,
                args.random_points,
                args.random_radius,
                args.seed,
                args.min_sep,
            )
            print(f"  random best:     lambda_min={lam:.12e} accepted={accepted} points={fmt_points(points)}")
        else:
            print("  random best:     skipped")
        print(f"  q-cache={len(qfun.cache)} elapsed={time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
