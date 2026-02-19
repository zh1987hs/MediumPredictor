from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .models import AnnotationRecord


def parse_provided_ko_table(path: str | Path) -> list[AnnotationRecord]:
    """Parse simple KO table: gene_id<TAB>ko (comma-separated allowed)."""
    records: list[AnnotationRecord] = []
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = row.get("gene_id") or row.get("gene") or row.get("id")
            ko_raw = row.get("ko") or row.get("KO") or ""
            if not gene_id:
                continue
            kos = [x.strip() for x in ko_raw.replace(";", ",").split(",") if x.strip()]
            records.append(AnnotationRecord(gene_id=gene_id, kos=kos, source="provided_ko", confidence=0.95))
    return records


def parse_generic_annotation(path: str | Path, source: str) -> list[AnnotationRecord]:
    """Best-effort parser for eggnog/dram-like TSV with gene_id and KO-like columns."""
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    gid_col = next((c for c in df.columns if c.lower() in {"gene_id", "query", "id", "locus_tag", "gene"}), df.columns[0])
    ko_col = next((c for c in df.columns if "ko" in c.lower() or "kegg_ko" in c.lower()), None)
    ec_col = next((c for c in df.columns if c.lower() in {"ec", "ec_number", "enzyme"}), None)

    records: list[AnnotationRecord] = []
    for _, row in df.iterrows():
        kos = []
        ecs = []
        if ko_col:
            kos = [x.strip().replace("ko:", "") for x in str(row[ko_col]).replace(";", ",").split(",") if x.strip()]
        if ec_col:
            ecs = [x.strip() for x in str(row[ec_col]).replace(";", ",").split(",") if x.strip()]
        records.append(AnnotationRecord(gene_id=str(row[gid_col]), kos=kos, ecs=ecs, source=source, confidence=0.8))
    return records


def run_prodigal(genome_fna: str | Path, proteins_out: str | Path) -> None:
    """Run prodigal to predict proteins from genome fasta."""
    cmd = ["prodigal", "-i", str(genome_fna), "-a", str(proteins_out), "-p", "single", "-q"]
    subprocess.run(cmd, check=True)


def parse_faa_ids(proteins_faa: str | Path) -> list[str]:
    """Extract protein IDs from FASTA as placeholder when annotations are missing."""
    return [rec.id for rec in SeqIO.parse(str(proteins_faa), "fasta")]
