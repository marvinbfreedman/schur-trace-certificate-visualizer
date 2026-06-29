#!/usr/bin/env python3
r"""Exact Mellin split identity for Xi atoms and Volterra atoms.

For one finite-core atom

    f_i(v) = a_i exp(beta_i v - c_i exp(v)),      v >= 0,

the positive-side Xi contribution is

    X_i(z) = 1/2 int_0^inf f_i(v) exp(i z v/2) dv
           = 1/2 a_i c_i^{-alpha} Gamma(alpha,c_i),
      alpha = beta_i + i z/2.

Splitting at the moving Volterra base point s gives the exact identity

    X_i(z) = B_i(s,z) + T_i(s,z),

where

    B_i(s,z) = 1/2 a_i c_i^{-alpha} int_{c_i}^{c_i e^s}
               y^{alpha-1} exp(-y) dy,

and

    T_i(s,z) = 1/2 Psi(s) exp(i z s/2)
               int_0^inf V_i(s,u) exp(i z u/2) du.

Here

    V_i(s,u) = rho_i(s) exp(beta_i u - c_i e^s(exp(u)-1))

is the corresponding normalized Volterra ratio atom inside
Psi(s+u)/Psi(s).

Thus the exact Mellin convolution does not give a constant 6x6 atom mixing
matrix.  It gives a diagonal moving-tail identity plus an s-dependent
incomplete-gamma boundary prefix.  This script verifies the identity and
measures the boundary/tail sizes on the de Branges sample nodes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / ".pydeps"))

import mpmath as mp  # noqa: E402

from klm_debranges_branch_transport_theorem import default_z_nodes  # noqa: E402
from klm_debranges_trace_map_constructor import to_float  # noqa: E402
from quotient_factorization_mp import endpoint_ratios, psi_value  # noqa: E402
from tilde3_volterra_confluent_mp import pieces_for  # noqa: E402
from xi_mellin_volterra_mode_match import ExactFiniteCoreXi, atom_records  # noqa: E402


def alpha(record: dict, z, fourier_scale: mp.mpf) -> mp.mpc:
    return mp.mpf(record["beta"]) + mp.j * fourier_scale * mp.mpc(z)


def atom_total(record: dict, z, fourier_scale: mp.mpf) -> mp.mpc:
    a = alpha(record, z, fourier_scale)
    c = mp.mpf(record["c"])
    return mp.mpf("0.5") * record["coeff"] * c ** (-a) * mp.gammainc(a, c, mp.inf)


def atom_prefix(record: dict, s: mp.mpf, z, fourier_scale: mp.mpf) -> mp.mpc:
    a = alpha(record, z, fourier_scale)
    c = mp.mpf(record["c"])
    return mp.mpf("0.5") * record["coeff"] * c ** (-a) * mp.gammainc(a, c, c * mp.e**s)


def atom_tail(record: dict, s: mp.mpf, z, fourier_scale: mp.mpf) -> mp.mpc:
    a = alpha(record, z, fourier_scale)
    c = mp.mpf(record["c"])
    return mp.mpf("0.5") * record["coeff"] * c ** (-a) * mp.gammainc(a, c * mp.e**s, mp.inf)


def volterra_atom_mellin(piece: tuple, ratio, s: mp.mpf, z, fourier_scale: mp.mpf) -> mp.mpc:
    _coeff, beta, c = piece
    a = mp.mpf(beta) + mp.j * fourier_scale * mp.mpc(z)
    lam = mp.mpf(c) * mp.e**s
    return ratio * mp.e**lam * lam ** (-a) * mp.gammainc(a, lam, mp.inf)


def identity_rows(args, xi_core: ExactFiniteCoreXi):
    pieces = pieces_for(args.kind)
    records = atom_records(args.kind)
    if len(pieces) != len(records):
        raise ValueError("piece/record mismatch")
    fourier_scale = mp.mpf(args.fourier_scale)
    omega = mp.mpf(args.omega)
    s_values = [mp.mpf(x) for x in args.s_values.replace(",", " ").split()]
    z_values = []
    for z in default_z_nodes():
        z_values.append(("base", z))
        z_values.append(("plus", mp.mpc(z) + mp.j * omega))
        z_values.append(("minus", mp.mpc(z) - mp.j * omega))

    rows = []
    max_total_error = mp.mpf("0")
    max_tail_error = mp.mpf("0")
    max_even_error = mp.mpf("0")
    max_boundary_fraction = mp.mpf("0")
    max_tail_fraction = mp.mpf("0")
    min_tail_fraction = mp.inf
    max_effective_multiplier = mp.mpf("0")
    for branch, z in z_values:
        even_total = mp.fsum(atom_total(record, z, fourier_scale) for record in records)
        even_core = xi_core(z) if args.fourier_scale == mp.mpf("0.5") else None
        if even_core is not None:
            max_even_error = max(
                max_even_error,
                abs(even_total - even_core) / max(abs(even_core), mp.mpf("1e-300")),
            )
        for s in s_values:
            psi_s = psi_value(s, pieces)
            ratios = endpoint_ratios(s, pieces)
            for idx, record in enumerate(records):
                total = atom_total(record, z, fourier_scale)
                prefix = atom_prefix(record, s, z, fourier_scale)
                tail = atom_tail(record, s, z, fourier_scale)
                mellin = volterra_atom_mellin(pieces[idx], ratios[idx][0], s, z, fourier_scale)
                tail_from_volterra = (
                    mp.mpf("0.5") * psi_s * mp.e ** (mp.j * fourier_scale * mp.mpc(z) * s) * mellin
                )
                total_error = abs(total - prefix - tail) / max(abs(total), mp.mpf("1e-300"))
                tail_error = abs(tail - tail_from_volterra) / max(abs(tail), mp.mpf("1e-300"))
                boundary_fraction = abs(prefix) / max(abs(total), mp.mpf("1e-300"))
                tail_fraction = abs(tail) / max(abs(total), mp.mpf("1e-300"))
                effective = (
                    2
                    * mp.e ** (-mp.j * fourier_scale * mp.mpc(z) * s)
                    * (total - prefix)
                    / max(abs(psi_s), mp.mpf("1e-300"))
                )
                effective_multiplier = abs(effective / total) if abs(total) else mp.inf
                max_total_error = max(max_total_error, total_error)
                max_tail_error = max(max_tail_error, tail_error)
                max_boundary_fraction = max(max_boundary_fraction, boundary_fraction)
                max_tail_fraction = max(max_tail_fraction, tail_fraction)
                min_tail_fraction = min(min_tail_fraction, tail_fraction)
                if effective_multiplier != mp.inf:
                    max_effective_multiplier = max(max_effective_multiplier, effective_multiplier)
                rows.append(
                    {
                        "branch": branch,
                        "z": {"real": to_float(mp.re(z)), "imag": to_float(mp.im(z))},
                        "s": to_float(s),
                        "atom": idx,
                        "mode": record["mode"],
                        "term": record["term"],
                        "totalAbs": to_float(abs(total)),
                        "prefixAbs": to_float(abs(prefix)),
                        "tailAbs": to_float(abs(tail)),
                        "boundaryFraction": to_float(boundary_fraction),
                        "tailFraction": to_float(tail_fraction),
                        "totalSplitRelativeError": to_float(total_error),
                        "tailVolterraRelativeError": to_float(tail_error),
                        "effectiveTailMultiplierAbs": to_float(effective_multiplier)
                        if effective_multiplier != mp.inf
                        else None,
                    }
                )
    return {
        "rows": rows,
        "maxTotalSplitRelativeError": to_float(max_total_error),
        "maxTailVolterraRelativeError": to_float(max_tail_error),
        "maxEvenAtomSumRelativeError": to_float(max_even_error),
        "maxBoundaryFraction": to_float(max_boundary_fraction),
        "minTailFraction": to_float(min_tail_fraction),
        "maxTailFraction": to_float(max_tail_fraction),
        "maxEffectiveTailMultiplierAbs": to_float(max_effective_multiplier),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["raw1", "raw2", "raw3", "tilde3"], default="tilde3")
    parser.add_argument("--omega", default="0.49")
    parser.add_argument("--fourier-scale", default="0.5")
    parser.add_argument("--s-values", default="0,0.02,0.1,0.3,0.545,1,2,4,8")
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--json-out", default="xi_mellin_convolution_boundary_identity.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    xi_core = ExactFiniteCoreXi(args.kind, "positive")
    result = identity_rows(args, xi_core)
    data = {
        "theoremName": "exact Mellin convolution boundary identity for finite Xi atoms",
        "kind": args.kind,
        "omega": float(args.omega),
        "fourierScale": float(args.fourier_scale),
        "identity": (
            "X_i(z)=B_i(s,z)+1/2 Psi(s) exp(i scale z s) "
            "int_0^inf V_i(s,u) exp(i scale z u)du"
        ),
        "conclusion": (
            "The exact atom transport is not a constant 6x6 mixing matrix. "
            "It is diagonal in the Volterra tail after multiplying by Psi(s), "
            "plus an s-dependent incomplete-gamma boundary prefix."
        ),
        **result,
        "nextProofTarget": (
            "Build the Hardy-to-Volterra transmutation with this moving boundary "
            "prefix: represent or cancel B_i(s,z) by endpoint/trace terms, then "
            "use the diagonal Volterra tail identity for the remaining tail."
        ),
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Exact Mellin convolution boundary identity")
    print(f"  kind={args.kind} omega={args.omega} scale={args.fourier_scale}")
    print(f"  max total split error: {data['maxTotalSplitRelativeError']:.3e}")
    print(f"  max tail/Volterra error: {data['maxTailVolterraRelativeError']:.3e}")
    print(f"  max atom-sum error: {data['maxEvenAtomSumRelativeError']:.3e}")
    print(f"  boundary fraction max: {data['maxBoundaryFraction']:.3e}")
    print(f"  tail fraction range: {data['minTailFraction']:.3e} .. {data['maxTailFraction']:.3e}")
    print(f"  max effective tail multiplier: {data['maxEffectiveTailMultiplierAbs']:.3e}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
