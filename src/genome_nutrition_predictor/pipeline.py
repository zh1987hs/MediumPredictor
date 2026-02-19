from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .annotation import parse_faa_ids, parse_generic_annotation, parse_provided_ko_table, run_prodigal
from .models import AnnotationRecord
from .reporting import write_outputs
from .rules import load_rules
from .scoring import index_ko_to_genes, infer_auxotrophy, score_capability, score_pathway


def load_annotations(annotation_mode: str, annotations: str | None, proteins: str | None, genome: str | None) -> list[AnnotationRecord]:
    """Load annotation records according to selected mode."""
    if annotation_mode == "provided_ko":
        if not annotations:
            raise ValueError("--annotations is required for provided_ko mode")
        return parse_provided_ko_table(annotations)

    if annotation_mode in {"eggnog_tsv", "dram"}:
        if not annotations:
            raise ValueError("--annotations required")
        return parse_generic_annotation(annotations, annotation_mode)

    if annotation_mode in {"prokka_gff", "bakta"}:
        if not annotations:
            raise ValueError("--annotations required")
        return parse_generic_annotation(annotations, annotation_mode)

    if annotation_mode == "kofam":
        raise NotImplementedError("kofam mode requires external KOfamScan integration; use provided_ko currently")

    if proteins:
        return [AnnotationRecord(gene_id=pid, kos=[], source="proteins_only", confidence=0.1) for pid in parse_faa_ids(proteins)]

    if genome:
        tmp = Path("predicted_proteins.faa")
        run_prodigal(genome, tmp)
        return [AnnotationRecord(gene_id=pid, kos=[], source="prodigal_only", confidence=0.1) for pid in parse_faa_ids(tmp)]

    raise ValueError("No supported inputs were provided")


def run_prediction(
    out_dir: str | Path,
    annotation_mode: str,
    annotations: str | None,
    proteins: str | None,
    genome: str | None,
    rules_path: str | None,
    db_dir: str | None,
    offline: bool,
    min_completeness_threshold: float,
    threads: int,
) -> dict[str, Any]:
    """Run full prediction pipeline for one genome/sample."""
    rules = load_rules(rules_path)
    ann = load_annotations(annotation_mode, annotations, proteins, genome)
    ko_index = index_ko_to_genes(ann)

    pathway_results = [score_pathway(p, ko_index) for p in rules.get("pathways", [])]
    by_id = {p.pathway_id: p for p in pathway_results}

    auxotrophy_rows: list[dict[str, Any]] = []
    for fac in rules.get("auxotrophy_factors", []):
        de_novo = by_id.get(fac["de_novo_pathway_id"])
        salvage = by_id.get(fac["salvage_pathway_id"]) if fac.get("salvage_pathway_id") else None
        if not de_novo:
            continue
        auxotrophy_rows.append(
            {
                "factor": fac["id"],
                "prediction": infer_auxotrophy(de_novo, salvage),
                "de_novo_status": de_novo.status,
                "salvage_status": salvage.status if salvage else None,
            }
        )

    cfg = rules.get("preference_model", {"transport_weight": 0.5, "catabolism_weight": 0.4, "regulation_weight": 0.1})
    capabilities = [score_capability(r, ko_index, cfg) for r in rules.get("substrate_rules", [])]

    metals = []
    for mr in rules.get("metal_rules", []):
        hit_kos = [ko for ko in mr.get("indicative_kos", []) if ko in ko_index]
        metals.append(
            {
                "factor": mr["id"],
                "message": mr["message_if_present"] if hit_kos else mr.get("message_if_absent", "No specific evidence"),
                "evidence_kos": hit_kos,
            }
        )

    metadata = {
        "tool": "genome-nutrition-predictor",
        "version": "0.1.0",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "annotation_mode": annotation_mode,
        "rules_path": rules_path or "rules/default_rules.yaml",
        "db_dir": db_dir,
        "offline": offline,
        "min_completeness_threshold": min_completeness_threshold,
        "threads": threads,
        "num_genes": len(ann),
    }

    write_outputs(out_dir, metadata, pathway_results, auxotrophy_rows, capabilities, metals)
    return {
        "metadata": metadata,
        "auxotrophy": auxotrophy_rows,
        "capabilities": [asdict(c) for c in capabilities],
    }
