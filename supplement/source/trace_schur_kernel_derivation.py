#!/usr/bin/env python3
r"""Derive the continuum trace-side Schur kernel.

This is a proof ledger, with one finite consistency check.

Let V be the completed Volterra/form domain, R:V->X_R the completed trace map,
and N=ker R.  Choose any section J:X_R->V with R J=I.  Write

    a(n,m)=Q(n,m),       n,m in N,
    b(n,x)=Q(n,Jx),
    c(x,y)=Q(Jx,Jy).

The closed quotient theorem gives a nonnegative closed form A for a and a
Douglas representative Gamma such that

    b(n,x)=<A^(1/2)n, Gamma x>.

The intrinsic trace-side Schur kernel is then

    D_trace(x,y)=c(x,y)-<Gamma x, Gamma y>.

This file records the exact derivation:

1. D_trace is independent of the section J.
2. D_trace is the relaxed constrained-energy kernel:

       D_trace(x,x) = inf { Q(f) : Rf=x }.

3. Therefore proving D_trace>=0 is equivalent to full-form positivity.

The last item is still the open analytic theorem.  The finite probes support it
but do not prove the continuum Volterra/Green Gram factorization.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp

from quotient_factorization_mp import (
    columns,
    diag,
    frob_norm,
    gram_matrix,
    positive_part_inverse,
    trace_matrix,
)


def load_json(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    out = {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
    }
    if blocker:
        out["blocker"] = blocker
    return out


def build_small_args(args: argparse.Namespace):
    return SimpleNamespace(
        model="full",
        kind="tilde3",
        omega=args.omega,
        L=args.L,
        basis=args.check_basis,
        quad=max(args.min_quad, args.quad_factor * args.check_basis),
        laguerre=args.laguerre,
        endpoint_kernel_order=args.endpoint_kernel_order,
        endpoint_kernel_rmax=args.endpoint_kernel_rmax,
        constraints=args.check_constraints,
        constraint_min=args.constraint_min,
        constraint_max=args.constraint_max,
        jet_order=args.jet_order,
        endpoint_order=args.endpoint_order,
        endpoint_rmax=args.endpoint_rmax,
        endpoint_tol=args.endpoint_tol,
        rank_tol=args.rank_tol,
        psd_tol=args.psd_tol,
        margin="0",
        dps=args.dps,
    )


def quotient_decomp(K: mp.matrix, R: mp.matrix, args):
    gram = R.T * R
    rvals, rvecs = mp.eigsy((gram + gram.T) / 2, eigvals_only=False)
    rmax = max(abs(v) for v in rvals) if len(rvals) else mp.mpf("0")
    rank_tol = mp.mpf(args.rank_tol) * max(mp.mpf("1"), rmax)
    n_idx = [i for i, val in enumerate(rvals) if val <= rank_tol]
    u_idx = [i for i, val in enumerate(rvals) if val > rank_tol]
    N = columns(rvecs, n_idx)
    U = columns(rvecs, u_idx)
    A = N.T * K * N
    B = N.T * K * U
    C = U.T * K * U
    _, _, _, _, Aplus = positive_part_inverse(A, mp.mpf(args.psd_tol))
    S = C - B.T * Aplus * B

    sigmas = [mp.sqrt(rvals[i]) for i in u_idx]
    RU = R * U
    W = mp.matrix(R.rows, U.cols)
    for j, sigma in enumerate(sigmas):
        for i in range(R.rows):
            W[i, j] = RU[i, j] / sigma
    sigma_inv = diag([1 / sigma for sigma in sigmas]) if sigmas else mp.matrix(0)
    D_trace = W * sigma_inv * S * sigma_inv * W.T if sigmas else mp.matrix(R.rows)
    return N, U, A, Aplus, B, C, (S + S.T) / 2, (D_trace + D_trace.T) / 2


def deterministic_section_change(ncols: int, ucols: int) -> mp.matrix:
    T = mp.matrix(ncols, ucols)
    for i in range(ncols):
        for j in range(ucols):
            sign = -1 if (i + 2 * j) % 2 else 1
            T[i, j] = sign * mp.mpf(i + 1) * mp.mpf(j + 2) / mp.mpf("37")
    return T


def finite_section_invariance_check(args: argparse.Namespace) -> dict:
    local = build_small_args(args)
    omega = mp.mpf(local.omega)
    length = mp.mpf(local.L)
    K, polys = gram_matrix(local, omega, length)
    _, R = trace_matrix(polys, local)
    N, U, A, Aplus, B, C, S, D_trace = quotient_decomp(K, R, local)

    T = deterministic_section_change(N.cols, U.cols)
    U2 = U + N * T
    B2 = N.T * K * U2
    C2 = U2.T * K * U2
    S2 = (C2 - B2.T * Aplus * B2 + (C2 - B2.T * Aplus * B2).T) / 2

    section_error = frob_norm(S2 - S)

    # Check transported trace quadratic identity:
    # y=RU z, y^T D_trace y = z^T S z.
    max_quad_error = mp.mpf("0")
    for seed in range(1, 5):
        z = mp.matrix(U.cols, 1)
        for i in range(U.cols):
            z[i, 0] = mp.mpf((seed + 1) * (i + 2)) / mp.mpf("11")
            if (seed + i) % 2:
                z[i, 0] *= -1
        y = R * U * z
        lhs = (y.T * D_trace * y)[0]
        rhs = (z.T * S * z)[0]
        max_quad_error = max(max_quad_error, abs(lhs - rhs))

    svals = mp.eigsy(S, eigvals_only=True)
    return {
        "basis": local.basis,
        "constraints": local.constraints,
        "traceRank": U.cols,
        "traceNullity": N.cols,
        "sectionChangeFrobeniusError": float(section_error),
        "traceTransportQuadraticMaxError": float(max_quad_error),
        "schurComplementMin": float(svals[0]) if len(svals) else 0.0,
        "schurComplementMax": float(svals[-1]) if len(svals) else 0.0,
        "closed": section_error < mp.mpf(args.check_tol)
        and max_quad_error < mp.mpf(args.check_tol),
    }


def finite_rows(*probes: dict) -> list[dict]:
    rows = []
    for probe in probes:
        for row in probe.get("rows", []):
            rows.append(
                {
                    "basis": row.get("basis"),
                    "constraints": row.get("constraints"),
                    "traceRank": row.get("traceRank"),
                    "schurComplementMin": row.get("schurComplementMin"),
                    "traceKernelMin": row.get("traceKernelMin"),
                    "weakModeDecay": (row.get("minModeExponentialFit") or {}).get("decay"),
                    "repairFreeFiniteClosed": row.get("repairFreeFiniteClosed"),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe-json", default="repair_free_schur_kernel_probe.json")
    parser.add_argument("--probe-b12-json", default="repair_free_schur_kernel_probe_b12c9.json")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--L", default="8")
    parser.add_argument("--check-basis", type=int, default=6)
    parser.add_argument("--check-constraints", type=int, default=3)
    parser.add_argument("--quad-factor", type=int, default=2)
    parser.add_argument("--min-quad", type=int, default=14)
    parser.add_argument("--laguerre", type=int, default=24)
    parser.add_argument("--endpoint-kernel-order", type=int, default=18)
    parser.add_argument("--endpoint-kernel-rmax", default="10")
    parser.add_argument("--constraint-min", default="0.02")
    parser.add_argument("--constraint-max", default="0.545")
    parser.add_argument("--jet-order", type=int, default=9)
    parser.add_argument("--endpoint-order", type=int, default=28)
    parser.add_argument("--endpoint-rmax", default="10")
    parser.add_argument("--endpoint-tol", default="1e-18")
    parser.add_argument("--rank-tol", default="1e-26")
    parser.add_argument("--psd-tol", default="1e-26")
    parser.add_argument("--check-tol", default="1e-35")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="trace_schur_kernel_derivation.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    probe = load_json(args.probe_json)
    probe_b12 = load_json(args.probe_b12_json)
    finite_check = finite_section_invariance_check(args)
    finite_evidence = finite_rows(probe, probe_b12)
    finite_closed = bool(finite_evidence) and all(
        bool(row.get("repairFreeFiniteClosed")) for row in finite_evidence
    )

    statuses = {
        "kernelDefinitionStatus": status(
            "intrinsic trace-side Schur kernel definition",
            True,
            (
                "For any section J of R, define "
                "D_trace(x,y)=Q(Jx,Jy)-<Gamma_J x,Gamma_J y>."
            ),
        ),
        "sectionIndependenceStatus": status(
            "section independence",
            finite_check["closed"],
            (
                "If J'=J+h with h:X_R->N, then Gamma' = Gamma + A^(1/2)h "
                "and the terms in c'-Gamma'^*Gamma' cancel algebraically. "
                "The finite section-change check also closes below tolerance."
            ),
            blocker=None if finite_check["closed"] else "Finite section-change check failed.",
        ),
        "infimumFormulaStatus": status(
            "constrained-energy infimum formula",
            True,
            (
                "Completing the square gives D_trace(x,x)=inf{Q(f): Rf=x} "
                "under the already imported Douglas range condition."
            ),
        ),
        "finiteTraceKernelStatus": status(
            "finite repair-free Schur kernels positive",
            finite_closed,
            (
                "All imported finite probes have min eig(C-Gamma^*Gamma)>0, "
                "including the 12/9 stress point."
            ),
            blocker=None if finite_closed else "A finite probe has negative Schur complement.",
        ),
        "volterraGreenGramStatus": status(
            "continuum Volterra/Green Gram positivity",
            False,
            (
                "The exact intrinsic kernel has been derived, but no explicit "
                "continuum Volterra/Green feature map has yet been constructed."
            ),
            blocker=(
                "Derive functions G_x(u) and a positive measure/operator dmu(u) "
                "such that D_trace(x,y)=integral G_x(u) G_y(u) dmu(u) in the "
                "transported X_R norm."
            ),
        ),
    }

    data = {
        "theoremName": "continuum trace-side Schur kernel derivation",
        "definitions": {
            "V": "completed Weyl/Volterra form domain",
            "R": "completed trace map V -> X_R",
            "N": "ker R",
            "J": "arbitrary section X_R -> V with R J = I",
            "A": "closed nonnegative form operator for Q restricted to N",
            "GammaJ": "Douglas representative b(n,x)=<A^(1/2)n,Gamma_J x>",
            "DTrace": "D_trace(x,y)=Q(Jx,Jy)-<Gamma_J x,Gamma_J y>",
        },
        "exactIdentities": [
            "D_trace is independent of the section J.",
            "D_trace(x,x)=inf{Q(f): Rf=x}.",
            "D_trace>=0 on X_R iff Q>=0 on the completed form domain.",
            (
                "In finite sampled trace coordinates RU=W Sigma, "
                "D_sample=W Sigma^{-1}(C-Gamma^*Gamma)Sigma^{-1}W^*."
            ),
        ],
        "statuses": statuses,
        "finiteSectionInvarianceCheck": finite_check,
        "finiteEvidenceRows": finite_evidence,
        "continuumTheoremClosed": False,
        "correctedConclusion": (
            "The exact continuum trace-side Schur kernel is now identified. "
            "The algebraic derivation and section invariance are closed. "
            "The positivity theorem remains the explicit Volterra/Green Gram "
            "factorization of this intrinsic kernel in the transported X_R norm."
        ),
        "nextProofTarget": (
            "Construct an explicit Green feature map for D_trace, probably by "
            "solving the Euler-Lagrange minimizer for inf{Q(f):Rf=x} and "
            "expressing the minimized energy as an integral of squares."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("trace-side Schur kernel derivation")
    print(f"  section invariance check: {finite_check['closed']}")
    print(f"  section error: {finite_check['sectionChangeFrobeniusError']:.3e}")
    print(f"  trace transport error: {finite_check['traceTransportQuadraticMaxError']:.3e}")
    print(f"  finite Schur kernels positive: {finite_closed}")
    print("  continuum Volterra/Green Gram positivity: open")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
