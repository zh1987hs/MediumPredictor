from __future__ import annotations

import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from .models import AnnotationRecord


def _split_flexible(line: str) -> list[str]:
    """Split a row by tab/comma/whitespace while keeping simple tokens."""
    return [x.strip().strip('"').strip("'") for x in re.split(r"[\t, ]+", line.strip()) if x.strip()]


def _norm_col(col: str) -> str:
    return col.strip().lstrip("\ufeff").lower()


def parse_provided_ko_table(path: str | Path) -> list[AnnotationRecord]:
    """Parse user KO table robustly.

    Supported formats:
    - tab-separated header: gene_id\tko
    - whitespace/comma-separated header/data: gene_id ko
    - files with UTF-8 BOM in header
    - headerless 2-column-like rows (gene_id KO)
    - duplicated rows collapsed by (gene_id, KO)
    """
    raw_lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8-sig", errors="ignore").splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not raw_lines:
        return []

    header_parts = _split_flexible(raw_lines[0])
    header_norm = [_norm_col(c) for c in header_parts]

    gid_idx: int | None = None
    ko_idx: int | None = None
    data_start = 1

    for i, c in enumerate(header_norm):
        if c in {"gene_id", "gene", "id", "query", "locus_tag"}:
            gid_idx = i
            break
    for i, c in enumerate(header_norm):
        if c in {"ko", "kegg_ko"} or "ko" in c:
            ko_idx = i
            break

    # If header not recognized, treat first line as data row and use first two columns.
    if gid_idx is None or ko_idx is None:
        gid_idx, ko_idx = 0, 1
        data_start = 0

    gene_to_kos: dict[str, set[str]] = defaultdict(set)
    for ln in raw_lines[data_start:]:
        parts = _split_flexible(ln)
        if max(gid_idx, ko_idx) >= len(parts):
            continue
        gene_id = parts[gid_idx].strip()
        ko_raw = parts[ko_idx].strip()
        if not gene_id or not ko_raw:
            continue

        # Support KO lists in one field (K00001;K00002 or comma-delimited)
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
