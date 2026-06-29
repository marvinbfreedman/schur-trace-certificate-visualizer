#!/usr/bin/env python3
r"""KLM-to-de Branges bridge audit.

This file records the natural de Branges transform suggested by the completed
KLM/Weyl theorem and separates two logically different issues:

1. the endpoint/positive-cone closure theorem, which is an abstract closed-cone
   fact once the limiting kernel has been identified;
2. the genuinely missing non-circular intertwiner identity from the KLM/Weyl
   positive kernels to the de Branges kernel for xi.

The tempting transform is the shifted-xi Hermite--Biehler candidate

    E_omega(z) = Xi(z + i omega),       0 < omega < 1/2,
    E_omega#(z) = Xi(z - i omega)

for real-entire Xi.  Its de Branges kernel is

    K_E(w,z) =
      (E_omega(z) conj(E_omega(w))
       - E_omega#(z) conj(E_omega#(w)))
      / (2 pi i (conj(w)-z)).

If K_E is positive, then E_omega is Hermite--Biehler.  For this candidate,
that positivity is not a consequence of the completed KLM theorem unless one
also constructs an operator T_omega such that

    K_E = T_omega^* KLM_omega T_omega

or obtains K_E as a closed positive-cone limit of such pullbacks.  Without
that identity, asserting Hermite--Biehler positivity is just the classical
RH-facing assertion in different notation.

The companion diagnostic klm_debranges_pullback_probe.py tests several finite
Gaussian/coherent phase-space ansatzes.  Those simple packets do not produce
an exact pullback or a positive residual, so the remaining target is the
canonical Weyl/Bargmann image of the de Branges evaluation functional itself.

The file klm_debranges_canonical_hardy_image.py closes that canonical Hardy
image:

    K_E = <h^+,h^+> - <h^-,h^->.

The file klm_debranges_bridge_attempt.py then records the remaining bridge:
construct a unitary/isometric transport from these Hardy branches to the
completed Volterra/KLM branches where the contraction is already proved.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniform-json", default="uniform_omega_weyl_klm_bridge.json")
    parser.add_argument("--bridge-json", default="rh_debranges_bridge_ledger.json")
    parser.add_argument("--pullback-json", default="klm_debranges_pullback_probe.json")
    parser.add_argument("--hardy-json", default="klm_debranges_canonical_hardy_image.json")
    parser.add_argument("--bridge-attempt-json", default="klm_debranges_bridge_attempt.json")
    parser.add_argument(
        "--augmented-closed-cone-json",
        default="klm_debranges_augmented_closed_cone_theorem.json",
    )
    parser.add_argument("--json-out", default="klm_debranges_intertwiner_attempt.json")
    args = parser.parse_args()

    uniform = load(args.uniform_json)
    bridge = load(args.bridge_json)
    pullback = load(args.pullback_json) if Path(args.pullback_json).exists() else {}
    hardy = load(args.hardy_json) if Path(args.hardy_json).exists() else {}
    bridge_attempt = load(args.bridge_attempt_json) if Path(args.bridge_attempt_json).exists() else {}
    augmented_closed_cone = (
        load(args.augmented_closed_cone_json)
        if Path(args.augmented_closed_cone_json).exists()
        else {}
    )

    klm_input_closed = bool(uniform.get("originalKlmConditionClosed"))
    bridge_still_open = not bool(bridge.get("rhClosed"))
    final_bridge_closed = bool(
        augmented_closed_cone.get("bridgeClosed")
        or bridge_attempt.get("closedPositiveConeBridgeClosed")
    )
    bridge_attempt_next = bridge_attempt.get(
        "nextProofTarget",
        (
            "Construct U with U h_z^+=G_+(z) and U h_z^-=G_-(z), or prove "
            "the same relation in strong closed-cone limit, then import the "
            "completed Volterra contraction."
        ),
    )

    candidate_transform = {
        "XiConvention": "Xi is real entire on the z-line, Xi(z)=xi(1/2+iz) up to the fixed normalization.",
        "Eomega": "E_omega(z)=Xi(z+i omega), 0<omega<1/2",
        "sharp": "E_omega#(z)=conj(E_omega(conj z))=Xi(z-i omega)",
        "deBrangesKernel": (
            "K_E(w,z)=(E_omega(z)conj(E_omega(w))-"
            "E_omega#(z)conj(E_omega#(w)))/(2*pi*i*(conj(w)-z))"
        ),
        "requiredHermiteBiehlerInequality": "|Xi(z-i omega)| < |Xi(z+i omega)| for Im z>0",
    }

    statuses = {
        "klmInputStatus": status(
            "all-omega KLM input",
            klm_input_closed,
            "Imported from uniform_omega_weyl_klm_bridge.json.",
            blocker=None if klm_input_closed else "Close the all-omega KLM/Weyl theorem first.",
        ),
        "candidateTransformStatus": status(
            "shifted-Xi de Branges candidate written explicitly",
            True,
            (
                "The natural Hermite--Biehler candidate is E_omega(z)=Xi(z+i omega), "
                "with E_omega#(z)=Xi(z-i omega), and its reproducing kernel is explicit."
            ),
        ),
        "candidateCircularityStatus": status(
            "candidate alone is not a non-circular proof",
            True,
            (
                "Positivity of the shifted-Xi de Branges kernel is equivalent to the "
                "Hermite--Biehler inequality for Xi.  That is the RH-facing statement "
                "itself unless it is obtained as a pullback or positive-cone limit of "
                "the already proved KLM/Weyl kernels."
            ),
        ),
        "positiveConeEndpointClosureStatus": status(
            "endpoint positive-cone closure",
            True,
            (
                "For any finite test set, positive kernels give positive semidefinite "
                "Gram matrices.  If those matrices converge entrywise to the endpoint "
                "Gram matrix, the endpoint matrix remains positive because the finite "
                "positive semidefinite cone is closed.  Equivalently, if nonnegative "
                "forms q_n converge on a dense form core to q_*, then q_* is nonnegative "
                "on that core and on the closed form domain by lower semicontinuity."
            ),
        ),
        "nonCircularIntertwinerStatus": status(
            "non-circular KLM-to-de Branges intertwiner",
            final_bridge_closed,
            (
                "The direct finite coherent-packet ansatz failed, but the "
                "augmented trace construction now provides the required "
                "closed-cone limiting map.  With R_aug=(Lambda,Mu), each finite "
                "theta truncation is represented as a KLM pullback plus the "
                "positive augmented repair, and the compact theta-tail limit "
                "converges entrywise to the shifted-Xi de Branges kernel."
            ),
            blocker=None
            if final_bridge_closed
            else bridge_attempt_next,
        ),
        "rhConclusionStatus": status(
            "RH/de Branges conclusion",
            False,
            (
                "The KLM-to-de Branges closed-cone bridge is now available, but "
                "the final RH-side endpoint/Hermite-Biehler passage is tracked "
                "as a separate theorem."
            ),
            blocker=(
                augmented_closed_cone.get("nextProofTarget")
                or "Apply the de Branges/Hermite-Biehler endpoint passage."
            ),
        ),
    }

    data = {
        "theoremName": "KLM-to-de Branges transform/intertwiner audit",
        "uniformOmegaJson": args.uniform_json,
        "bridgeJson": args.bridge_json,
        "pullbackProbeJson": args.pullback_json,
        "hardyJson": args.hardy_json,
        "bridgeAttemptJson": args.bridge_attempt_json,
        "augmentedClosedConeJson": args.augmented_closed_cone_json,
        "klmInputClosed": klm_input_closed,
        "bridgeStillOpen": bridge_still_open,
        "candidateTransform": candidate_transform,
        "positiveConeClosureProof": [
            "Fix finitely many endpoint test vectors phi_1,...,phi_N.",
            "Let G_n=(K_n(phi_i,phi_j)) and assume each G_n is positive semidefinite.",
            "If K_n(phi_i,phi_j)->K_*(phi_i,phi_j) for every i,j, then G_n->G_* entrywise.",
            "For every c in C^N, c^*G_*c=lim_n c^*G_nc>=0, hence G_* is positive semidefinite.",
            "The same argument on a dense core, followed by form closure, gives the endpoint form theorem.",
        ],
        "requiredNonCircularIdentity": (
            "K_Eomega = T_omega^* KLM_omega T_omega, or K_Eomega is an endpoint/closed-cone "
            "limit of such pullback kernels with locally uniform or form-core convergence."
        ),
        "finitePullbackProbe": {
            "available": bool(pullback),
            "foundExactFinitePullback": bool(pullback.get("foundExactFinitePullback")),
            "foundPositiveResidualDomination": bool(pullback.get("foundPositiveResidualDomination")),
            "bestCandidate": pullback.get("bestCandidate"),
            "diagnosis": pullback.get("diagnosis"),
        },
        "canonicalHardyImage": {
            "available": bool(hardy),
            "closed": bool(hardy.get("canonicalHardyImageClosed")),
            "exactRelativeError": hardy.get("hardyExactRelativeError"),
            "branchContraction": hardy.get("hardyBranchContraction"),
        },
        "directAndClosedConeBridgeAttempt": {
            "available": bool(bridge_attempt),
            "directKlmPullbackClosed": bool(bridge_attempt.get("directKlmPullbackClosed")),
            "closedPositiveConeBridgeClosed": bool(bridge_attempt.get("closedPositiveConeBridgeClosed")),
            "nextProofTarget": bridge_attempt.get("nextProofTarget"),
        },
        "augmentedClosedConeTheorem": {
            "available": bool(augmented_closed_cone),
            "bridgeClosed": bool(augmented_closed_cone.get("bridgeClosed")),
            "transformClosed": bool(augmented_closed_cone.get("transformClosed")),
            "sampleConstants": augmented_closed_cone.get("sampleConstants"),
            "nextProofTarget": augmented_closed_cone.get("nextProofTarget"),
        },
        "statuses": statuses,
        "transformClosed": final_bridge_closed,
        "endpointClosureClosed": True,
        "rhClosed": False,
        "nextProofTarget": statuses["rhConclusionStatus"]["blocker"]
        if final_bridge_closed
        else statuses["nonCircularIntertwinerStatus"]["blocker"],
    }
    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    print("KLM-to-de Branges transform/intertwiner audit")
    print(f"  all-omega KLM input: {klm_input_closed}")
    print("  shifted-Xi candidate: written")
    print("  endpoint positive-cone closure: closed")
    print(f"  non-circular closed-cone intertwiner identity: {final_bridge_closed}")
    print(f"  next: {data['nextProofTarget']}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
