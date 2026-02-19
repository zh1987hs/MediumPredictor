from pathlib import Path

from genome_nutrition_predictor.annotation import parse_provided_ko_table


def test_parse_provided_ko_whitespace_and_duplicates(tmp_path: Path):
    p = tmp_path / "input.tsv"
    p.write_text(
        "gene_id ko\n"
        "cxk-1000001 K02313\n"
        "cxk-1000001 K02313\n"
        "cxk-1000002 K02338\n"
        "cxk-1000002 K02338\n"
        "cxk-1000002 K02338\n"
    )
    rows = parse_provided_ko_table(p)
    by_gene = {r.gene_id: r.kos for r in rows}
    assert by_gene["cxk-1000001"] == ["K02313"]
    assert by_gene["cxk-1000002"] == ["K02338"]
    assert len(rows) == 2
