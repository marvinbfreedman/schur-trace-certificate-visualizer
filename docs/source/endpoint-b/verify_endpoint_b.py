#!/usr/bin/env python3
"""Certify the fixed negative direction for the defining endpoint-B kernel."""

from flint import arb, ctx


# Exact decimals obtained by rounding the archived witness to 50 decimal places.
DIRECTION = (
    "-0.00068032524193180803154064513095808308689565764525",
    "0.00839120211421558324858907779775057389972563512764",
    "-0.04825918681658664312617000928564889447072005844862",
    "0.17135113839037381677686404931369720758091899926691",
    "-0.41905941422734026350039921191792416703251770258178",
    "0.74538910758315986188825453494822102270354146357036",
    "-0.99199063908159270541221361123216815853547177510003",
    "1",
    "-0.76349718019770744209903196210720488589152822174732",
    "0.43598695093804101765268313087939822801665565353143",
    "-0.18096510110124778634972104830104893185541012650673",
    "0.05169279645538153491126079330057912755145273701376",
    "-0.00910637758991081800301235656611612652489444808522",
    "0.00074702877506803986783166465980564646393453861360",
)


def moment(p, a):
    """Enclose the integral of x**p exp(-a*(x-1)) over x >= 1."""
    return (a - (p + 1) * a.log()).exp() * a.gamma_upper(p + 1)


def log_moment(p, a, step):
    center = (moment(p + step, a) - moment(p - step, a)) / (2 * step)
    error = step**2 * moment(p + step + 1, a)
    return center + arb(0, error.abs_upper())


def product_moment(alpha, beta, moment_at):
    # Expand (1-alpha*x)(1-beta*x) against x**(3/2).
    a = alpha + beta
    p = arb(3) / 2
    return moment_at(p, a) - a * moment_at(p + 1, a) + alpha * beta * moment_at(p + 2, a)


def main():
    ctx.prec = 768
    step = arb(2) ** -120
    base = arb.pi()
    nodes = [arb(3) / 5 + index * arb(3) / 50 for index in range(len(DIRECTION))]
    coeffs = [arb(value) for value in DIRECTION]
    if not (arb("0.545") < nodes[0] and nodes[-1] < 8):
        raise ArithmeticError("The witness nodes do not avoid the active trace interval")

    rates = [base, *(base * node.exp() for node in nodes)]
    ordinary = {}
    logarithmic = {}
    for i, alpha in enumerate(rates):
        for j, beta in enumerate(rates[: i + 1]):
            ordinary[i, j] = ordinary[j, i] = product_moment(alpha, beta, moment)
            logarithmic[i, j] = logarithmic[j, i] = product_moment(
                alpha, beta, lambda p, a: log_moment(p, a, step)
            )

    quadratic = arb(0)
    for i, s in enumerate(nodes, 1):
        for j, t in enumerate(nodes, 1):
            log_part = (
                logarithmic[i, j] - logarithmic[i, 0]
                - logarithmic[0, j] + logarithmic[0, 0]
            )
            cross_part = (
                t * (ordinary[i, j] - ordinary[0, j])
                + s * (ordinary[i, j] - ordinary[i, 0])
            ) / 2
            entry = (log_part + cross_part) / ((rates[i] - base) * (rates[j] - base))
            quadratic += coeffs[i - 1] * coeffs[j - 1] * entry

    if not (arb("-1.119e-28") < quadratic and quadratic < arb("-1.118e-28")):
        raise ArithmeticError(f"Expected a strictly negative enclosure, obtained {quadratic}")
    print(f"endpoint-B quadratic form at {ctx.prec} bits: {quadratic}")
    print("PASS: witness is exterior to the active trace interval and strictly negative")


if __name__ == "__main__":
    main()
