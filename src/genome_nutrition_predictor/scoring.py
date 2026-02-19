from __future__ import annotations

from collections import defaultdict
from typing import Any

from .models import AnnotationRecord, CapabilityResult, PathwayResult


def index_ko_to_genes(records: list[AnnotationRecord]) -> dict[str, list[str]]:
    """Build KO -> supporting gene IDs index."""
    idx: dict[str, list[str]] = defaultdict(list)
    for rec in records:
        for ko in rec.kos:
            idx[ko].append(rec.gene_id)
    return idx


def _step_hit(step: dict[str, Any], ko_index: dict[str, list[str]]) -> tuple[bool, list[str], list[str]]:
    candidates = step.get("any_of", [])
    hit_kos = [ko for ko in candidates if ko in ko_index]
    genes: list[str] = []
    for ko in hit_kos:
        genes.extend(ko_index.get(ko, []))
    return bool(hit_kos), sorted(set(hit_kos)), sorted(set(genes))


def score_pathway(pathway: dict[str, Any], ko_index: dict[str, list[str]]) -> PathwayResult:
    """Score a rule-based pathway with essential/optional weighted completeness."""
    essential_steps = pathway.get("essential_steps", [])
    optional_steps = pathway.get("optional_steps", [])

    essential_hits = 0
    optional_hits = 0
    missing_ess: list[str] = []
    evidence_genes: set[str] = set()

    for step in essential_steps:
        hit, _, genes = _step_hit(step, ko_index)
        if hit:
            essential_hits += 1
            evidence_genes.update(genes)
        else:
            missing_ess.append(step.get("id", "unnamed"))

    for step in optional_steps:
        hit, _, genes = _step_hit(step, ko_index)
        if hit:
            optional_hits += 1
            evidence_genes.update(genes)

    essential_ratio = essential_hits / len(essential_steps) if essential_steps else 1.0
    optional_ratio = optional_hits / len(optional_steps) if optional_steps else 1.0

    ew = pathway.get("essential_weight", 0.8)
    ow = pathway.get("optional_weight", 0.2)
    completeness = ew * essential_ratio + ow * optional_ratio
    score = round(completeness * 100, 2)

    if essential_ratio >= 0.95 and completeness >= 0.8:
        status = "Present"
    elif essential_ratio >= 0.4 or completeness >= 0.4:
        status = "Partial"
    else:
        status = "Absent"

    return PathwayResult(
        pathway_id=pathway["id"],
        name=pathway.get("name", pathway["id"]),
        category=pathway.get("category", "unknown"),
        status=status,
        score_0_100=score,
        essential_ratio=round(essential_ratio, 3),
        optional_ratio=round(optional_ratio, 3),
        missing_essential_groups=missing_ess,
        evidence_genes=sorted(evidence_genes),
    )


def infer_auxotrophy(pathway: PathwayResult, salvage_pathway: PathwayResult | None = None) -> str:
    """Infer likely requirement label from de novo and salvage states."""
    salvage_present = salvage_pathway is not None and salvage_pathway.status in {"Present", "Partial"}
    if pathway.status == "Present":
        return "not required"
    if pathway.status == "Absent" and salvage_present:
        return "likely required (salvage evidence)"
    if pathway.status == "Absent":
        return "likely required"
    return "possible required"


def score_capability(rule: dict[str, Any], ko_index: dict[str, list[str]], cfg: dict[str, float]) -> CapabilityResult:
    """Evaluate substrate capability and preference from transport + catabolism rules."""
    transport_steps = rule.get("transport_steps", [])
    catabolic_steps = rule.get("catabolic_steps", [])
    regulation_steps = rule.get("regulation_markers", [])

    def ratio(steps: list[dict[str, Any]]) -> tuple[float, set[str], set[str]]:
        if not steps:
            return 1.0, set(), set()
        hits = 0
        genes: set[str] = set()
        kos: set[str] = set()
        for s in steps:
            hit, hit_kos, hit_genes = _step_hit(s, ko_index)
            if hit:
                hits += 1
                genes.update(hit_genes)
                kos.update(hit_kos)
        return hits / len(steps), genes, kos

    tr_ratio, tr_genes, tr_kos = ratio(transport_steps)
    ca_ratio, ca_genes, ca_kos = ratio(catabolic_steps)
    rg_ratio, rg_genes, rg_kos = ratio(regulation_steps)

    capability = (tr_ratio >= 0.5 or rule.get("passive_diffusion", False)) and ca_ratio >= 0.6
    pref = cfg["transport_weight"] * tr_ratio + cfg["catabolism_weight"] * ca_ratio + cfg["regulation_weight"] * rg_ratio

    genes = tr_genes | ca_genes | rg_genes
    kos = tr_kos | ca_kos | rg_kos
    return CapabilityResult(
        substrate=rule["id"],
        kind=rule.get("kind", "carbon"),
        capability=capability,
        preference_score=round(min(pref, 1.0), 3),
        transport_score=round(tr_ratio, 3),
        catabolism_score=round(ca_ratio, 3),
        regulation_score=round(rg_ratio, 3),
        evidence_gene_ids=sorted(genes),
        evidence_kos=sorted(kos),
    )
