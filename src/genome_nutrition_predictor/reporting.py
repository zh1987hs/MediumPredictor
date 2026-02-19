from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pandas as pd

from .models import CapabilityResult, PathwayResult


def write_outputs(
    out_dir: str | Path,
    metadata: dict[str, Any],
    pathway_results: list[PathwayResult],
    auxotrophy: list[dict[str, Any]],
    capabilities: list[CapabilityResult],
    metals: list[dict[str, Any]],
) -> None:
    """Write required output files for one genome."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    (out / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False))

    pathways_df = pd.DataFrame([asdict(x) for x in pathway_results])
    pathways_df.to_csv(out / "pathways.tsv", sep="\t", index=False)

    evidence_rows = []
    for pr in pathway_results:
        for gid in pr.evidence_genes:
            evidence_rows.append({"context": pr.pathway_id, "gene_id": gid, "status": pr.status})
    for cp in capabilities:
        for gid in cp.evidence_gene_ids:
            evidence_rows.append({"context": cp.substrate, "gene_id": gid, "status": "capability" if cp.capability else "weak"})
    pd.DataFrame(evidence_rows).to_csv(out / "evidence.tsv", sep="\t", index=False)

    medium_df = build_medium_suggestion(auxotrophy, capabilities)
    medium_df.to_csv(out / "medium_suggestion.tsv", sep="\t", index=False)

    summary = {
        "metadata": metadata,
        "auxotrophy": auxotrophy,
        "capabilities": [asdict(x) for x in capabilities],
        "metals": metals,
        "pathway_counts": pathways_df["status"].value_counts().to_dict() if not pathways_df.empty else {},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    report = render_report(summary, pathway_results)
    (out / "report.md").write_text(report)


def build_medium_suggestion(auxotrophy: list[dict[str, Any]], capabilities: list[CapabilityResult]) -> pd.DataFrame:
    """Build compact minimal-medium recommendation table."""
    required = [a["factor"] for a in auxotrophy if a["prediction"].startswith("likely required")]
    possible = [a["factor"] for a in auxotrophy if a["prediction"].startswith("possible required")]

    carbon_top = sorted([c for c in capabilities if c.kind == "carbon" and c.capability], key=lambda x: x.preference_score, reverse=True)[:3]
    nitrogen_top = sorted([c for c in capabilities if c.kind == "nitrogen" and c.capability], key=lambda x: x.preference_score, reverse=True)[:3]

    return pd.DataFrame(
        [
            {"section": "base", "item": "M9-like salts + trace elements", "reason": "default basal medium"},
            {"section": "supplement_required", "item": ", ".join(required) or "none", "reason": "de novo absent"},
            {"section": "supplement_possible", "item": ", ".join(possible) or "none", "reason": "de novo partial"},
            {"section": "carbon_top3", "item": ", ".join([c.substrate for c in carbon_top]) or "unknown", "reason": "transport+catabolism score"},
            {"section": "nitrogen_top3", "item": ", ".join([n.substrate for n in nitrogen_top]) or "unknown", "reason": "transport+catabolism score"},
        ]
    )


def render_report(summary: dict[str, Any], pathway_results: list[PathwayResult]) -> str:
    """Create markdown report."""
    lines = ["# Genome Nutrition Predictor Report", "", "## Auxotrophy / Requirement"]
    for row in summary["auxotrophy"]:
        lines.append(f"- **{row['factor']}**: {row['prediction']} (de novo={row['de_novo_status']}, salvage={row.get('salvage_status', 'NA')})")

    lines.extend(["", "## Substrate capability and preference"])
    for cap in sorted(summary["capabilities"], key=lambda x: x["preference_score"], reverse=True):
        lines.append(f"- {cap['kind']}::{cap['substrate']}: capability={cap['capability']}, score={cap['preference_score']}")

    lines.extend(["", "## Metals and trace element evidence"])
    for metal in summary["metals"]:
        lines.append(f"- {metal['factor']}: {metal['message']}")

    lines.extend(["", "## Pathway overview"])
    for p in pathway_results:
        lines.append(f"- {p.category}::{p.name} => {p.status} ({p.score_0_100})")
    lines.append("")
    return "\n".join(lines)
