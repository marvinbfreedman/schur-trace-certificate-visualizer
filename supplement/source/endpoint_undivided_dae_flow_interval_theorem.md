# Undivided DAE Graph and Endpoint-Flow Collocation Interval Theorem

## Statement

On the closed interval `[0.02, 0.545]`, let the bounded seven-plane be written
as `X=T(q)Y`, where `T(q)=(I_7;q)`. For the piecewise polynomial coefficient
field named in `endpoint_undivided_dae_krawczyk_enclosure.json`, the equations

```text
a8(q' + qB + (qb)q) + c7 q + d = 0,
G' = (B+bq)G
```

have unique nodal solutions in the displayed Arb collocation boxes. The graph
equation is never divided by `a8`; its crossing collocation chart therefore
remains defined at the certified simple zero. The induced collocation endpoint matrix

```text
M = C_R T(q_R)G_R - C_L T(q_L)G_L
```

has active row rank two in this polynomial collocation model.

## Verification

For each graph chart, a high-precision Newton center is enclosed by

```text
K(X) = x0 - R F(x0) + (I - R F'(X))(X-x0).
```

All seven graph boxes and all seven induced-flow boxes are outward-rounded Arb
Krawczyk self-maps. Boundary uncertainty is transported from one flow chart to
the next. The best endpoint determinant uses columns `(0,1)` and satisfies

```text
|det M_(0,1)| >= 250055428.246444.
```

The complete radii and contractions are recorded in the two dependency JSON
files and are recomputed by their corresponding Python scripts.

## Scope

This proves the finite polynomial collocation theorem.  The separate exact
crossing theorem now supplies a uniform between-node Bernstein/Cauchy
enclosure for `e^(q)(s)`, `q<=8`, and a continuous residual-to-solution
inclusion on the singular crossing interval.  Full continuum endpoint
promotion still requires the corresponding exact regular-piece propagation
and transport of its error tube into the restricted active endpoint minor.
The theorem JSON therefore records
`exactCrossingContinuumPromotionClosed = true` while correctly retaining
`exactContinuumTheoremClosed = false`.
