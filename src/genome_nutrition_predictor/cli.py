from __future__ import annotations

import json
from pathlib import Path

import typer

from .pipeline import run_prediction
from .rules import build_template, load_rules

app = typer.Typer(help="Predict likely nutrient requirements and substrate preferences from bacterial genome annotations")


@app.command()
def run(
    out: str = typer.Option(..., help="Output directory"),
    genome: str | None = typer.Option(None, help="Genome .fna for CDS prediction"),
    proteins: str | None = typer.Option(None, help="Protein FASTA"),
    annotations: str | None = typer.Option(None, help="Annotation table"),
    annotation_mode: str = typer.Option("provided_ko", help="kofam/provided_ko/eggnog_tsv/dram/prokka_gff/bakta"),
    rules: str | None = typer.Option(None, help="Custom rules YAML"),
    db_dir: str | None = typer.Option(None, help="Database directory"),
    offline: bool = typer.Option(False, help="Disable online download attempts"),
    min_completeness_threshold: float = typer.Option(0.8, help="Completeness threshold"),
    threads: int = typer.Option(1, help="Threads"),
    log_level: str = typer.Option("INFO", help="Log level"),
) -> None:
    """Run prediction for a single genome/sample."""
    _ = log_level
    result = run_prediction(
        out_dir=out,
        annotation_mode=annotation_mode,
        annotations=annotations,
        proteins=proteins,
        genome=genome,
        rules_path=rules,
        db_dir=db_dir,
        offline=offline,
        min_completeness_threshold=min_completeness_threshold,
        threads=threads,
    )
    typer.echo(json.dumps(result["metadata"], indent=2, ensure_ascii=False))


@app.command("build-rules-template")
def build_rules_template(out: str = typer.Option(..., help="Output YAML path")) -> None:
    """Write a starter YAML for custom rule extension."""
    Path(out).write_text(build_template())
    typer.echo(f"Template written to: {out}")


@app.command()
def explain(factor: str, rules: str | None = typer.Option(None)) -> None:
    """Explain pathway rules for a nutrient factor/substrate."""
    rs = load_rules(rules)
    fac = [x for x in rs.get("auxotrophy_factors", []) if x["id"].lower() == factor.lower()]
    sub = [x for x in rs.get("substrate_rules", []) if x["id"].lower() == factor.lower()]
    if fac:
        typer.echo(json.dumps(fac[0], indent=2, ensure_ascii=False))
        return
    if sub:
        typer.echo(json.dumps(sub[0], indent=2, ensure_ascii=False))
        return
    typer.echo(f"No rule matched factor={factor}")


@app.command()
def batch(
    input_list: str = typer.Option(..., help="TSV with columns: sample,annotations,proteins,genome,annotation_mode"),
    out: str = typer.Option(..., help="Output root"),
    rules: str | None = typer.Option(None),
    threads: int = typer.Option(1),
) -> None:
    """Batch processing of multiple genomes."""
    root = Path(out)
    root.mkdir(parents=True, exist_ok=True)
    lines = Path(input_list).read_text().strip().splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        sample = row["sample"]
        run_prediction(
            out_dir=root / sample,
            annotation_mode=row.get("annotation_mode", "provided_ko"),
            annotations=row.get("annotations") or None,
            proteins=row.get("proteins") or None,
            genome=row.get("genome") or None,
            rules_path=rules,
            db_dir=None,
            offline=True,
            min_completeness_threshold=0.8,
            threads=threads,
        )
    typer.echo(f"Batch complete: {out}")


if __name__ == "__main__":
    app()
