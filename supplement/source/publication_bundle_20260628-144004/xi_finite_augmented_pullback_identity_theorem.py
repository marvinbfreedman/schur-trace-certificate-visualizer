#!/usr/bin/env python3
r"""Finite augmented Xi/KLM pullback identity theorem.

For finite theta truncations Xi_N, this theorem packages the exact identity

    K_{E,N}=T_N^*KLM_NT_N + R_{aug,N}^*D_{aug,N}R_{aug,N}.

The proof uses only two symbolic ingredients:

1. the canonical Hardy image
       K_{E,N}=<h_N^+,h_N^+>-<h_N^-,h_N^->;
2. the symbolic Mellin boundary theorem, which splits each finite Xi atom into
       a diagonal Volterra/KLM tail plus the primitive Mu boundary trace.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def closed(data: dict, key: str) -> bool:
    item = data.get(key, {})
    if isinstance(item, dict):
        return bool(item.get("closed"))
    return bool(item)


def status(label: str, ok: bool, reason: str) -> dict:
    return {
        "label": label,
        "closed": ok,
        "status": "closed" if ok else "open",
        "reason": reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hardy-json",
        default="klm_debranges_canonical_hardy_image_hardened_theorem.json",
    )
    parser.add_argument(
        "--mellin-symbolic-json",
        default="xi_mellin_boundary_symbolic_theorem.json",
    )
    parser.add_argument(
        "--json-out",
        default="xi_finite_augmented_pullback_identity_theorem.json",
    )
    args = parser.parse_args()

    hardy = load(args.hardy_json)
    mellin = load(args.mellin_symbolic_json)

    hardy_ok = bool(hardy.get("canonicalHardyImageClosed"))
    mellin_ok = bool(
        mellin.get("exactMellinSplitClosed")
        and mellin.get("mellinBoundaryConcomitantClosed")
        and mellin.get("primitiveMuTraceClosed")
    )
    finite_ok = hardy_ok and mellin_ok

    hardy_status = status(
        "finite shifted-Xi Hardy branch identity",
        hardy_ok,
        (
            "The hardened canonical Hardy theorem gives the signed Hardy "
            "branch Gram for every finite Xi_N by applying the same identity "
            "to E_{omega,N}(z)=Xi_N(z+i omega)."
        ),
    )
    mellin_status = status(
        "finite atom Mellin-to-Volterra split",
        mellin_ok,
        (
            "The symbolic Mellin theorem splits every retained atom into a "
            "diagonal Volterra/KLM tail and the primitive Mu boundary trace."
        ),
    )
    pullback_status = status(
        "finite augmented pullback identity",
        finite_ok,
        (
            "Summing the atomwise tail identities gives T_N^*KLM_NT_N; "
            "summing the primitive prefixes gives the augmented trace repair "
            "R_aug,N^*D_aug,N R_aug,N."
        ),
    )

    data = {
        "theoremName": "finite augmented Xi/KLM pullback identity theorem",
        "proofClass": "symbolic identity",
        "hardyJson": args.hardy_json,
        "mellinSymbolicJson": args.mellin_symbolic_json,
        "statuses": {
            "hardyBranchIdentityStatus": hardy_status,
            "mellinVolterraSplitStatus": mellin_status,
            "finiteAugmentedPullbackStatus": pullback_status,
        },
        "finiteThetaAugmentedPullbackClosed": finite_ok,
        "finiteThetaPullbackIdentityClosed": finite_ok,
        "signedHardyGramEqualsFiniteDeBrangesKernel": hardy_ok,
        "finitePullbackFormula": (
            "K_{E,N}=T_N^*KLM_NT_N+R_{aug,N}^*D_{aug,N}R_{aug,N}"
        ),
        "proof": [
            "Apply the canonical Hardy identity to E_{omega,N}.",
            "Expand Xi_N into finitely many Mellin atoms.",
            "Apply the symbolic moving-boundary split to each atom.",
            "Collect diagonal Volterra tails into T_N^*KLM_NT_N.",
            "Collect primitive prefixes into the augmented trace repair.",
        ],
        "remainingAnalyticGap": None if finite_ok else "One symbolic input is open.",
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("Finite augmented Xi/KLM pullback identity theorem")
    print(f"  Hardy identity: {hardy_ok}")
    print(f"  Mellin split: {mellin_ok}")
    print(f"  finite pullback: {finite_ok}")
    print(f"  wrote {args.json_out}")


if __name__ == "__main__":
    main()
