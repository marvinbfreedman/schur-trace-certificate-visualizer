#!/usr/bin/env python3
r"""Build a publication-risk dependency DAG for the RH/KLM certificate ledgers.

This is an audit tool, not a proof checker.  It freezes the current formal
claim, follows JSON theorem/certificate imports, classifies dependency edges
conservatively, and ranks the links that need publication-grade review first.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, deque
from pathlib import Path
from typing import Any


CLAIM = (
    "KLM/Weyl positivity -> augmented closed-cone de Branges kernel positivity "
    "-> shifted-Xi endpoint zero exclusion"
)

DANGER_ZONES = {
    "original_weyl_kernel_positivity_assembly_theorem.json",
    "xi_augmented_trace_continuum_lift.json",
    "klm_debranges_augmented_closed_cone_theorem.json",
}

BREAK_QUESTIONS = [
    "Is the domain exactly the same on both sides?",
    "Are closures/density limits valid and explicitly proved?",
    "Is positivity strict or only semidefinite, and is strictness ever used?",
    "Are boundary terms actually zero in the completed domain?",
    "Is any finite certificate promoted to continuum without a compactness theorem?",
]

THEOREM_NAME_RE = re.compile(
    r"(theorem|certificate|audit|ledger|bridge|endpoint|augmented|uniform|"
    r"closure|continuum|klm|weyl|debranges|schur|proof)",
    re.IGNORECASE,
)

JSON_REF_RE = re.compile(r"[\w./-]+\.json\b")
TOKEN_RE = re.compile(r"[A-Z]?[a-z]+|[A-Z]+(?=[A-Z]|$)|\d+")
GENERIC_TOKENS = {
    "json",
    "xi",
    "klm",
    "weyl",
    "de",
    "branges",
    "augmented",
    "trace",
    "theorem",
    "certificate",
    "bridge",
    "status",
    "input",
    "closed",
    "closure",
}

EDGE_CLASSES = [
    "analytic proof",
    "symbolic identity",
    "interval/ball certificate",
    "numerical evidence",
    "unproven assumption",
]


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def as_text(value: Any, *, limit: int = 120_000) -> str:
    try:
        text = json.dumps(value, sort_keys=True, ensure_ascii=True)
    except TypeError:
        text = str(value)
    return text[:limit]


def is_theorem_like(path: Path, data: dict[str, Any]) -> bool:
    if "theoremName" in data or "statuses" in data:
        return True
    return bool(THEOREM_NAME_RE.search(path.name))


def is_closed_value(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, dict):
        if "closed" in value and isinstance(value["closed"], bool):
            return value["closed"]
        if value.get("status") == "closed":
            return True
        if value.get("status") == "open":
            return False
    return None


def collect_closed_flags(data: dict[str, Any]) -> dict[str, bool]:
    flags: dict[str, bool] = {}
    for key, value in data.items():
        if key == "statuses" and isinstance(value, dict):
            for skey, svalue in value.items():
                closed = is_closed_value(svalue)
                if closed is not None:
                    flags[f"statuses.{skey}"] = closed
            continue
        if key.endswith("Closed") or key in {
            "bridgeClosed",
            "transformClosed",
            "rhClosed",
            "formalRhClosed",
            "conditionalRhClosed",
            "endpointPassageClosed",
            "continuumAugmentedRepairClosed",
            "originalKlmConditionClosed",
            "originalWeylKernelPositivityClosed",
            "uniformOmegaCoverageClosed",
            "rhFacingChainClosed",
        }:
            closed = is_closed_value(value)
            if closed is not None:
                flags[key] = closed
    return flags


def recursively_find_json_refs(value: Any, path: tuple[str, ...] = ()) -> list[tuple[str, str]]:
    refs: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = (*path, str(key))
            if isinstance(child, str):
                for match in JSON_REF_RE.findall(child):
                    refs.append((".".join(child_path), Path(match).name))
            elif isinstance(child, list):
                for idx, item in enumerate(child):
                    refs.extend(recursively_find_json_refs(item, (*child_path, str(idx))))
            elif isinstance(child, dict):
                refs.extend(recursively_find_json_refs(child, child_path))
    elif isinstance(value, list):
        for idx, item in enumerate(value):
            refs.extend(recursively_find_json_refs(item, (*path, str(idx))))
    elif isinstance(value, str):
        for match in JSON_REF_RE.findall(value):
            refs.append((".".join(path), Path(match).name))
    return refs


def choose_context(source: dict[str, Any], target_name: str, import_key: str) -> str:
    chunks: list[Any] = [import_key, source.get("theoremName"), source.get("proof"), source.get("formalProof")]
    statuses = source.get("statuses")
    if isinstance(statuses, dict):
        target_stem = Path(target_name).stem
        import_leaf = import_key.split(".")[-1]
        tokens = {
            token.lower()
            for text in [target_stem, import_leaf]
            for token in TOKEN_RE.findall(text.replace("_", " "))
            if len(token) > 3 and token.lower() not in GENERIC_TOKENS
        }
        for key, value in statuses.items():
            key_l = key.lower()
            value_text = as_text(value, limit=10_000).lower()
            if any(token in key_l or token in value_text for token in tokens):
                chunks.append({key: value})
    for key in ("nextProofTarget", "caution", "scope", "notClaimedHere"):
        if key in source:
            chunks.append({key: source[key]})
    if "sample" in import_key.lower() or "pullback_limit" in target_name.lower():
        for key in ("sampleConstants", "bestCase"):
            if key in source:
                chunks.append({key: source[key]})
    return "\n".join(as_text(chunk, limit=20_000) for chunk in chunks if chunk is not None)


def classify_edge(import_key: str, source: dict[str, Any], target: dict[str, Any] | None, target_name: str) -> tuple[str, str, list[str]]:
    context = choose_context(source, target_name, import_key)
    if target:
        context += "\n" + target_name + "\n" + as_text(
            {
                "theoremName": target.get("theoremName"),
                "statuses": target.get("statuses"),
                "proof": target.get("proof"),
                "formalProof": target.get("formalProof"),
                "caution": target.get("caution"),
            },
            limit=40_000,
        )
    text = context.lower()
    reasons: list[str] = []

    unproven_hits = [
        "unproven",
        "not proved",
        "not claimed",
        "missing",
        "open",
        "blocker",
        "conditional only",
        "external proof",
        "independent external proof vetted\": false",
    ]
    numerical_hits = [
        "sample",
        "scan",
        "observed",
        "finite run",
        "relative error",
        "roundoff",
        "basis",
        "galerkin",
        "numerical",
        "probe",
        "smoke",
        "least-squares",
    ]
    interval_hits = [
        "interval",
        "ball",
        "bernstein",
        "krawczyk",
        "enclosure",
        "radius",
        "radii",
        "capacity",
        "quadrature error",
    ]
    symbolic_hits = [
        "identity",
        "exact",
        "symbolic",
        "mellin split",
        "green identity",
        "concomitant",
        "diagonal inequality",
        "normalization",
    ]
    analytic_hits = [
        "analytic",
        "proof",
        "theorem",
        "closure",
        "density",
        "compactness",
        "lower semicontinuity",
        "continuum",
        "bounded",
        "closed form",
        "mosco",
    ]

    if target is None:
        return "unproven assumption", "unknown", ["referenced JSON is missing or unreadable"]

    source_proof_class = str(source.get("proofClass", "")).lower()
    target_proof_class = str(target.get("proofClass", "")).lower()
    source_theorem_name = str(source.get("theoremName", "")).lower()
    edge_label = f"{import_key} {target_name}".lower()
    if (
        ("symbolic identity" in source_proof_class or "symbolic identity" in target_proof_class)
        and any(word in edge_label for word in ["consequence", "wrapper", "identity"])
    ):
        reasons.append("explicit symbolic consequence/wrapper proof metadata")
        return "symbolic identity", "certificate", reasons
    if (
        "symbolic identity" in source_proof_class
        and any(word in source_theorem_name for word in ["input", "consequence", "wrapper"])
    ):
        reasons.append("explicit symbolic input/consequence wrapper source metadata")
        return "symbolic identity", "certificate", reasons

    if any(hit in text for hit in unproven_hits) and not any(hit in text for hit in interval_hits + analytic_hits):
        reasons.append("open/missing/unproven language without a stronger certificate keyword")
        return "unproven assumption", "unknown", reasons

    if any(hit in text for hit in numerical_hits):
        if any(hit in text for hit in interval_hits):
            reasons.append("numerical language paired with interval/ball certificate terms")
            return "interval/ball certificate", "certificate", reasons
        if any(
            hit in text
            for hit in [
                "on sample",
                "sample verification",
                "closedconebridgeclosedonsample",
                "signedhardygramequalsdebrangeskernelonsample",
            ]
        ):
            reasons.append("explicit sample-only verification language")
            return "numerical evidence", "numerical", reasons
        if any(hit in text for hit in analytic_hits) and any(
            hit in text
            for hit in [
                "continuum theorem",
                "compactness theorem",
                "mosco",
                "lower semicontinuity",
                "weierstrass",
                "uniform compact convergence",
            ]
        ):
            reasons.append("finite/numerical evidence is paired with continuum closure language")
            return "analytic proof", "certificate", reasons
        reasons.append("sample/scan/Galerkin/error language without rigorous enclosure in edge context")
        return "numerical evidence", "numerical", reasons

    if any(hit in text for hit in interval_hits):
        reasons.append("interval/ball/enclosure terminology found")
        return "interval/ball certificate", "certificate", reasons
    if any(hit in text for hit in symbolic_hits):
        reasons.append("exact identity/symbolic/concomitant terminology found")
        return "symbolic identity", "certificate", reasons
    if any(hit in text for hit in analytic_hits):
        reasons.append("analytic proof/closure/continuum terminology found")
        return "analytic proof", "certificate", reasons

    reasons.append("no decisive proof metadata found")
    return "unproven assumption", "unknown", reasons


def risk_for_edge(edge_class: str, confidence: str, source: str, target: str, context: str) -> tuple[int, list[str]]:
    risk = 0
    reasons: list[str] = []
    if edge_class == "unproven assumption":
        risk += 100
        reasons.append("unproven assumption edge")
    elif edge_class == "numerical evidence":
        risk += 85
        reasons.append("numerical-only evidence")
    elif edge_class == "analytic proof":
        risk += 35
        reasons.append("analytic proof requires human audit")
    elif edge_class == "symbolic identity":
        risk += 25
        reasons.append("symbolic identity requires normalization audit")
    elif edge_class == "interval/ball certificate":
        risk += 20
        reasons.append("machine certificate requires radius/enclosure audit")

    if confidence == "unknown":
        risk += 20
        reasons.append("unknown confidence")
    if source in DANGER_ZONES or target in DANGER_ZONES:
        risk += 30
        reasons.append("danger-zone dependency")

    text = context.lower()
    symbolic_wrapper_edge = (
        edge_class == "symbolic identity"
        and (
            "explicit symbolic consequence/wrapper proof metadata" in text
            or "explicit symbolic input/consequence wrapper source metadata" in text
            or "consequence" in source.lower()
            or "consequence" in target.lower()
        )
    )
    if symbolic_wrapper_edge and source not in DANGER_ZONES and target not in DANGER_ZONES:
        reasons.append("symbolic consequence wrapper edge")
        return risk, reasons

    if any(word in text for word in ["domain", "space", "trace-fiber", "completed", "quotient"]):
        risk += 10
        reasons.append("domain/space matching issue")
    if any(word in text for word in ["closure", "density", "limit", "continuum", "galerkin"]):
        risk += 12
        reasons.append("closure/density/continuum passage issue")
    if "semidefinite" in text or "strict" in text:
        risk += 8
        reasons.append("strict vs semidefinite positivity issue")
    if any(word in text for word in ["boundary", "endpoint", "concomitant", "annihilat", "vanish"]):
        risk += 10
        reasons.append("boundary term issue")
    if any(word in text for word in ["finite", "sample", "scan", "observed", "basis"]):
        risk += 12
        reasons.append("finite-to-continuum promotion issue")
    return risk, reasons


def risk_for_node(name: str, data: dict[str, Any], outgoing_edges: list[dict[str, Any]]) -> tuple[int, list[str]]:
    risk = 0
    reasons: list[str] = []
    flags = collect_closed_flags(data)
    if name in DANGER_ZONES:
        risk += 40
        reasons.append("danger-zone node")
    if not flags:
        risk += 15
        reasons.append("no explicit closed flags")
    elif any(value is False for value in flags.values()):
        risk += 35
        reasons.append("contains open/false closed flags")

    text = as_text(data, limit=80_000).lower()
    if "independent external proof vetted\": false" in text:
        risk += 25
        reasons.append("explicit external-vetting caveat")
    if any(word in text for word in ["sample", "scan", "observed", "basis", "relative error"]):
        risk += 20
        reasons.append("contains numerical evidence language")
    if any(word in text for word in ["closure", "density", "continuum", "galerkin", "compactness"]):
        risk += 15
        reasons.append("contains closure/continuum passage language")
    if any(word in text for word in ["boundary", "endpoint", "concomitant", "annihilat", "vanish"]):
        risk += 12
        reasons.append("contains boundary/endpoint language")
    if outgoing_edges:
        risk += max(edge["riskScore"] for edge in outgoing_edges) // 3
    return risk, reasons


def reachable_from(root: str, edges_by_source: dict[str, list[dict[str, Any]]]) -> set[str]:
    seen: set[str] = set()
    queue: deque[str] = deque([root])
    while queue:
        current = queue.popleft()
        if current in seen:
            continue
        seen.add(current)
        for edge in edges_by_source.get(current, []):
            target = edge["target"]
            if target not in seen:
                queue.append(target)
    return seen


def dot_graph(nodes: dict[str, dict[str, Any]], edges: list[dict[str, Any]], reachable: set[str]) -> str:
    lines = ["digraph publication_audit {", "  rankdir=LR;"]
    for name in sorted(reachable):
        node = nodes.get(name)
        if node:
            label = f"{name}\\nrisk={node['riskScore']}"
            color = "red" if name in DANGER_ZONES else ("orange" if node["riskScore"] >= 70 else "gray")
        else:
            label = f"{name}\\nimported artifact"
            color = "lightgray"
        lines.append(f'  "{name}" [label="{label}", color="{color}"];')
    for edge in edges:
        if edge["source"] in reachable and edge["target"] in reachable:
            label = f"{edge['edgeClass']}\\nrisk={edge['riskScore']}"
            lines.append(f'  "{edge["source"]}" -> "{edge["target"]}" [label="{label}"];')
    lines.append("}")
    return "\n".join(lines)


def markdown_report(data: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# Publication Audit Dependency Graph")
    lines.append("")
    lines.append("## Frozen Claim")
    lines.append("")
    lines.append(data["frozenClaim"])
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- JSON files scanned: {data['counts']['jsonFilesScanned']}")
    lines.append(f"- Theorem/certificate nodes: {data['counts']['theoremLikeNodes']}")
    lines.append(f"- Reachable nodes from root: {data['counts']['reachableNodes']}")
    lines.append(f"- Edges: {data['counts']['edges']}")
    lines.append(f"- Formal chain closed in ledger: {data['formalChainClosed']}")
    lines.append(f"- Independent external proof vetted: {data['independentExternalProofVetted']}")
    lines.append("")
    lines.append("## Edge Classes")
    lines.append("")
    for key, value in data["counts"]["edgeClasses"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## Top Blockers")
    lines.append("")
    for item in data["rankedBlockers"][:15]:
        lines.append(
            f"- risk {item['riskScore']}: `{item['source']}` -> `{item['target']}` "
            f"({item['edgeClass']}, {item['confidence']})"
        )
        if item["riskReasons"]:
            lines.append(f"  - reasons: {', '.join(item['riskReasons'][:5])}")
        if item["classificationReasons"]:
            lines.append(f"  - classification: {', '.join(item['classificationReasons'][:3])}")
    lines.append("")
    lines.append("## Danger-Zone Audit")
    lines.append("")
    for name, audit in data["dangerZoneAudits"].items():
        lines.append(f"### `{name}`")
        lines.append("")
        lines.append(f"- risk: {audit['nodeRiskScore']}")
        lines.append(f"- closed flags: {audit['closedFlags']}")
        lines.append("- break questions:")
        for question in audit["breakQuestions"]:
            lines.append(f"  - {question}")
        if audit["topIncomingEdges"]:
            lines.append("- top incoming edges:")
            for edge in audit["topIncomingEdges"]:
                lines.append(
                    f"  - `{edge['source']}` -> `{edge['target']}` "
                    f"risk {edge['riskScore']} ({edge['edgeClass']})"
                )
        if audit["topOutgoingEdges"]:
            lines.append("- top outgoing edges:")
            for edge in audit["topOutgoingEdges"]:
                lines.append(
                    f"  - `{edge['source']}` -> `{edge['target']}` "
                    f"risk {edge['riskScore']} ({edge['edgeClass']})"
                )
        lines.append("")
    lines.append("## Downgrade Rules")
    lines.append("")
    for rule in data["downgradeRules"]:
        lines.append(f"- {rule}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="rh_formal_conclusion_consequence_theorem.json")
    parser.add_argument("--json-out", default="publication_audit_dependency_graph.json")
    parser.add_argument("--md-out", default="publication_audit_dependency_graph.md")
    args = parser.parse_args()

    cwd = Path.cwd()
    json_files = sorted(
        path
        for path in cwd.glob("*.json")
        if path.name != Path(args.json_out).name
        and not path.name.startswith("publication_audit_dependency_graph.")
    )
    raw: dict[str, dict[str, Any]] = {}
    for path in json_files:
        data = load_json(path)
        if data is not None:
            raw[path.name] = data

    theorem_like = {
        name: data for name, data in raw.items() if is_theorem_like(Path(name), data)
    }

    edges: list[dict[str, Any]] = []
    for source_name, source_data in theorem_like.items():
        seen_refs: set[tuple[str, str]] = set()
        for import_key, target_name in recursively_find_json_refs(source_data):
            if target_name == source_name or not target_name.endswith(".json"):
                continue
            if (import_key, target_name) in seen_refs:
                continue
            seen_refs.add((import_key, target_name))
            target_data = raw.get(target_name)
            edge_class, confidence, class_reasons = classify_edge(
                import_key, source_data, target_data, target_name
            )
            context = choose_context(source_data, target_name, import_key)
            risk_score, risk_reasons = risk_for_edge(
                edge_class, confidence, source_name, target_name, context
            )
            edges.append(
                {
                    "source": source_name,
                    "target": target_name,
                    "importKey": import_key,
                    "edgeClass": edge_class,
                    "confidence": confidence,
                    "riskScore": risk_score,
                    "riskReasons": risk_reasons,
                    "classificationReasons": class_reasons,
                    "targetExists": target_data is not None,
                }
            )

    edges_by_source: dict[str, list[dict[str, Any]]] = {}
    incoming_by_target: dict[str, list[dict[str, Any]]] = {}
    for edge in edges:
        edges_by_source.setdefault(edge["source"], []).append(edge)
        incoming_by_target.setdefault(edge["target"], []).append(edge)

    reachable = reachable_from(args.root, edges_by_source) if args.root in theorem_like else set()

    nodes: dict[str, dict[str, Any]] = {}
    for name, data in theorem_like.items():
        node_edges = edges_by_source.get(name, [])
        risk_score, risk_reasons = risk_for_node(name, data, node_edges)
        flags = collect_closed_flags(data)
        nodes[name] = {
            "path": name,
            "theoremName": data.get("theoremName"),
            "closedFlags": flags,
            "allExplicitClosedFlagsTrue": bool(flags) and all(flags.values()),
            "nextProofTarget": data.get("nextProofTarget"),
            "riskScore": risk_score,
            "riskReasons": risk_reasons,
            "importCount": len(node_edges),
            "importedJsons": sorted({edge["target"] for edge in node_edges}),
            "reachableFromRoot": name in reachable,
            "isDangerZone": name in DANGER_ZONES,
        }

    ranked_edges = sorted(
        [edge for edge in edges if edge["source"] in reachable],
        key=lambda item: (-item["riskScore"], item["source"], item["target"]),
    )

    danger_audits: dict[str, Any] = {}
    for name in sorted(DANGER_ZONES):
        node = nodes.get(name)
        danger_audits[name] = {
            "present": node is not None,
            "nodeRiskScore": node["riskScore"] if node else None,
            "nodeRiskReasons": node["riskReasons"] if node else ["missing danger-zone node"],
            "closedFlags": node["closedFlags"] if node else {},
            "breakQuestions": BREAK_QUESTIONS,
            "topIncomingEdges": sorted(
                incoming_by_target.get(name, []),
                key=lambda item: -item["riskScore"],
            )[:5],
            "topOutgoingEdges": sorted(
                edges_by_source.get(name, []),
                key=lambda item: -item["riskScore"],
            )[:5],
        }

    edge_class_counts = Counter(edge["edgeClass"] for edge in edges)
    root_data = raw.get(args.root, {})
    data = {
        "theoremName": "publication audit dependency graph",
        "frozenClaim": CLAIM,
        "root": args.root,
        "formalChainClosed": bool(root_data.get("formalRhClosed") or root_data.get("rhClosed")),
        "independentExternalProofVetted": bool(root_data.get("independentRhProofVetted")),
        "counts": {
            "jsonFilesScanned": len(json_files),
            "jsonFilesLoaded": len(raw),
            "theoremLikeNodes": len(theorem_like),
            "reachableNodes": len(reachable),
            "edges": len(edges),
            "reachableEdges": len(ranked_edges),
            "edgeClasses": {key: edge_class_counts.get(key, 0) for key in EDGE_CLASSES},
        },
        "nodes": nodes,
        "edges": edges,
        "reachableFromRoot": sorted(reachable),
        "rankedBlockers": ranked_edges[:50],
        "dangerZoneAudits": danger_audits,
        "downgradeRules": [
            "Sampling, scans, Galerkin basis evidence, or relative-error checks are numerical evidence unless paired with interval/ball or explicit continuum-closure language.",
            "Missing JSON imports or open/blocker language become unproven assumptions.",
            "Closure, density, quotient-domain, endpoint, and boundary-term claims always increase risk even when marked closed.",
            "Danger-zone nodes and their incident edges receive extra risk weight.",
            "The audit may report the formal chain closed while still ranking publication blockers.",
        ],
        "dotGraph": dot_graph(nodes, edges, reachable),
    }

    Path(args.json_out).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    Path(args.md_out).write_text(markdown_report(data), encoding="utf-8")

    print("Publication audit dependency graph")
    print(f"  frozen claim: {CLAIM}")
    print(f"  json files scanned: {len(json_files)}")
    print(f"  theorem/certificate nodes: {len(theorem_like)}")
    print(f"  reachable nodes from {args.root}: {len(reachable)}")
    print(f"  edges: {len(edges)}")
    print(f"  formal chain closed: {data['formalChainClosed']}")
    print(f"  independent external proof vetted: {data['independentExternalProofVetted']}")
    print("  edge classes:")
    for key in EDGE_CLASSES:
        print(f"    {key}: {edge_class_counts.get(key, 0)}")
    if ranked_edges:
        top = ranked_edges[0]
        print(
            "  top blocker: "
            f"{top['source']} -> {top['target']} "
            f"risk={top['riskScore']} class={top['edgeClass']}"
        )
    print(f"  wrote {args.json_out}")
    print(f"  wrote {args.md_out}")


if __name__ == "__main__":
    main()
