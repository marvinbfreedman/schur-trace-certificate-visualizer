#!/usr/bin/env python3
r"""Audit the primitive transport boundary terms.

For a half-line parity kernel P with mixed derivative H=P_xy and decay

    P(x,y) = int_x^inf int_y^inf H(u,v) du dv,

the original coordinate quadratic form is transported by the cumulative
primitive

    F(u) = int_0^u f(x) dx

as

    int int P(x,y) f(x) g(y) dx dy
      = int int H(u,v) F(u) G(v) du dv.

Equivalently, integrating H(u,v)F(u)G(v) by parts in u and v has no endpoint
contribution: F(0)=G(0)=0 kills the lower endpoints and the decay of P kills
the upper endpoints.  This means the primitive Weyl/Volterra transport itself
does not produce a canonical positive trace repair.  In the notation of the
boundary-comparison ledger,

    beta_bdy = Q_original - Q_Phi = 0

for the already identified mixed/Volterra form Q_Phi.

Therefore the previous sufficient comparison D_bdy >= D_q reduces to the much
sharper condition D_q=0 on the primitive trace image.  If that cannot be
proved, one must prove positivity of Q_Phi directly on the primitive image
rather than adding a nonexistent boundary repair.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_optional(path: str) -> dict | None:
    candidate = Path(path)
    if not candidate.exists():
        return None
    return json.loads(candidate.read_text(encoding="utf-8"))


def read(path: str) -> str:
    return Path(path).read_text(encoding="utf-8")


def status(label: str, closed: bool, reason: str, *, blocker: str | None = None) -> dict:
    return {
        "label": label,
        "closed": closed,
        "status": "closed" if closed else "open",
        "reason": reason,
        "blocker": blocker,
    }


def has_all(text: str, needles: list[str]) -> bool:
    return all(needle in text for needle in needles)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--notes", default="rh_klm_notes.md")
    parser.add_argument("--draft", default="rh_weyl_positive_draft.tex")
    parser.add_argument(
        "--quotient-json",
        default="weyl_volterra_quotient_schur_theorem.json",
    )
    parser.add_argument(
        "--boundary-comparison-json",
        default="canonical_boundary_repair_comparison.json",
    )
    parser.add_argument("--json-out", default="primitive_boundary_transport_audit.json")
    args = parser.parse_args()

    notes = read(args.notes)
    draft = read(args.draft)
    quotient = load_optional(args.quotient_json)
    comparison = load_optional(args.boundary_comparison_json)

    primitive_identity_recorded = has_all(
        notes,
        [
            "P_pm(x,y)=int_x^inf int_y^inf H_pm(u,v) du dv",
            "F(u)=int_0^u f(x) dx",
            "Volterra/log normalization then transports this mixed-kernel quadratic form",
        ],
    ) and has_all(
        draft,
        [
            "P_\\pm(x,y)=\\int_x^\\infty\\int_y^\\infty H_\\pm(u,v)",
            "F(u)=\\int_0^u f(x)\\,dx",
        ],
    )

    quotient_closed = bool(
        quotient and quotient.get("globalWeylVolterraSchurStatus", {}).get("closed")
    )
    dq_defined = bool(
        comparison
        and comparison.get("statuses", {})
        .get("minimalQuotientRepairStatus", {})
        .get("closed")
    )

    primitive_transport = status(
        "primitive transport has no endpoint boundary term",
        primitive_identity_recorded,
        (
            "For compact smooth half-line tests, F(0)=G(0)=0.  Since the parity "
            "kernels decay as either variable tends to infinity, integration "
            "by parts in the identity H=P_xy gives no lower or upper endpoint "
            "contribution.  Hence the coordinate form equals the mixed "
            "Volterra form exactly."
        ),
        blocker=None
        if primitive_identity_recorded
        else "Record the half-line primitive identity and the endpoint decay assumptions.",
    )
    beta_zero = status(
        "canonical primitive boundary form is zero",
        primitive_transport["closed"],
        (
            "With Q_Phi defined as the mixed/Volterra transport of the original "
            "coordinate form, beta_bdy=Q_original-Q_Phi is identically zero.  "
            "Thus D_bdy is the zero operator on the completed trace range."
        ),
        blocker=None
        if primitive_transport["closed"]
        else "First close the no-endpoint primitive transport identity.",
    )
    descends = status(
        "zero boundary form descends through R_global",
        beta_zero["closed"],
        (
            "The zero form annihilates ker R_global and is bounded in every "
            "transported trace norm, so it descends trivially with D_bdy=0."
        ),
    )
    comparison_reduction = status(
        "D_bdy >= D_q reduces to D_q=0",
        beta_zero["closed"] and dq_defined,
        (
            "Because D_bdy=0, the sufficient comparison D_bdy>=D_q is "
            "equivalent to D_q=0 on the relevant completed primitive trace "
            "range.  A positive hidden boundary repair cannot be supplied by "
            "the primitive Green identity."
        ),
        blocker=None
        if beta_zero["closed"] and dq_defined
        else "Define both D_bdy and D_q in the same trace range.",
    )
    dq_zero = status(
        "minimal quotient repair vanishes on the primitive trace image",
        False,
        (
            "This has not been proved.  It is now the exact replacement for "
            "the failed boundary-repair route.  One must either prove D_q R F=0 "
            "for every primitive image F, or prove Q_Phi>=0 directly on that "
            "primitive image."
        ),
        blocker=(
            "Characterize the primitive trace image Y=R({F: F'=original test}) "
            "inside X_R and prove D_q|_Y=0, or prove the equivalent direct "
            "positivity of Q_Phi on the primitive image."
        ),
    )
    original_positivity = status(
        "original positivity from primitive transport",
        False,
        (
            "The primitive transport is exact, but the quotient Schur theorem "
            "only gives automatic positivity on ker R_global.  Since primitive "
            "images are not contained in ker R_global and D_bdy=0, original "
            "positivity still requires D_q to vanish on the primitive trace "
            "image or a direct primitive-image positivity theorem."
        ),
        blocker=dq_zero["blocker"],
    )

    data = {
        "theoremName": "primitive boundary transport audit",
        "notes": args.notes,
        "draft": args.draft,
        "quotientJson": args.quotient_json if quotient else None,
        "boundaryComparisonJson": args.boundary_comparison_json if comparison else None,
        "importedSignals": {
            "primitiveIdentityRecorded": primitive_identity_recorded,
            "quotientSchurClosed": quotient_closed,
            "minimalDqDefined": dq_defined,
        },
        "statuses": {
            "primitiveTransportNoEndpointStatus": primitive_transport,
            "canonicalPrimitiveBoundaryZeroStatus": beta_zero,
            "zeroBoundaryDescendsStatus": descends,
            "comparisonReducesToDqZeroStatus": comparison_reduction,
            "dqZeroOnPrimitiveImageStatus": dq_zero,
            "originalPositivityFromPrimitiveTransportStatus": original_positivity,
        },
        "integrationByPartsProof": [
            "Start from P(x,y)=int_x^inf int_y^inf H(u,v)dudv, so H=P_xy.",
            "Let F(u)=int_0^u f(x)dx and G(v)=int_0^v g(y)dy.",
            "Then int int HFG = [P_yFG]_{u=0}^{inf}-int int P_y f G.",
            "The boundary term vanishes because F(0)=0 and P_y(inf,v)=0.",
            "A second integration by parts gives int int P f g; the v-boundary vanishes because G(0)=0 and P(x,inf)=0.",
        ],
        "operators": {
            "Dbdy": "0 on X_R for the primitive transport already identified with Q_Phi",
            "Dq": "minimal quotient/Douglas repair from the closed Schur theorem",
        },
        "correctedNextProofTarget": dq_zero["blocker"],
        "primitiveBoundaryTransportAuditClosed": (
            primitive_transport["closed"] and beta_zero["closed"] and descends["closed"]
        ),
        "originalWeylKernelPositivityClosed": False,
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Primitive boundary transport audit")
    print(f"  primitive identity recorded: {primitive_identity_recorded}")
    print(f"  no endpoint boundary term: {primitive_transport['closed']}")
    print(f"  D_bdy=0: {beta_zero['closed']}")
    print(f"  D_bdy>=D_q reduces to D_q=0: {comparison_reduction['closed']}")
    print(f"  D_q zero on primitive image: {dq_zero['closed']}")
    print(f"  next: {data['correctedNextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
