from __future__ import annotations

import csv
import subprocess
from collections import defaultdict
from pathlib import Path

from .models import AnnotationRecord


def parse_provided_ko_table(path: str | Path) -> list[AnnotationRecord]:
    """Parse user KO table robustly.

    Supported formats:
    - tab-separated header: gene_id\tko
    - whitespace-separated header/data: gene_id ko
    - duplicated rows are collapsed by (gene_id, KO)
    """
    raw_lines = [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not raw_lines:
        return []

    header = raw_lines[0]
    data_lines = raw_lines[1:]

    # Decide parser mode: explicit TSV first; fallback to whitespace split.
    use_tsv = "\t" in header

    gene_to_kos: dict[str, set[str]] = defaultdict(set)

    if use_tsv:
        reader = csv.DictReader(raw_lines, delimiter="\t")
        for row in reader:
            gene_id = (row.get("gene_id") or row.get("gene") or row.get("id") or "").strip()
            ko_raw = (row.get("ko") or row.get("KO") or "").strip()
            if not gene_id:
                continue
            for ko in [x.strip() for x in ko_raw.replace(";", ",").split(",") if x.strip()]:
                gene_to_kos[gene_id].add(ko)
    else:
        cols = header.split()
        if len(cols) < 2:
            return []
        try:
            gid_idx = next(i for i, c in enumerate(cols) if c.lower() in {"gene_id", "gene", "id", "query", "locus_tag"})
            ko_idx = next(i for i, c in enumerate(cols) if c.lower() in {"ko", "kegg_ko"} or "ko" in c.lower())
        except StopIteration:
            return []

        for ln in data_lines:
            parts = ln.split()
            if max(gid_idx, ko_idx) >= len(parts):
                continue
            gene_id = parts[gid_idx].strip()
            ko_raw = parts[ko_idx].strip()
            if not gene_id:
                continue
            for ko in [x.strip() for x in ko_raw.replace(";", ",").split(",") if x.strip()]:
                gene_to_kos[gene_id].add(ko)

    return [
        AnnotationRecord(gene_id=gid, kos=sorted(kos), source="provided_ko", confidence=0.95)
        for gid, kos in sorted(gene_to_kos.items())
    ]


def parse_generic_annotation(path: str | Path, source: str) -> list[AnnotationRecord]:
    """Best-effort parser for eggnog/dram-like TSV with gene_id and KO-like columns."""
    import pandas as pd
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
    from Bio import SeqIO
    return [rec.id for rec in SeqIO.parse(str(proteins_faa), "fasta")]
