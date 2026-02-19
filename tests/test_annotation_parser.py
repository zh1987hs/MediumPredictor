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


def test_parse_provided_ko_with_bom_header(tmp_path: Path):
    p = tmp_path / "bom.tsv"
    p.write_text("\ufeffgene_id ko\nabc K00001\n")
    rows = parse_provided_ko_table(p)
    assert len(rows) == 1
    assert rows[0].gene_id == "abc"
    assert rows[0].kos == ["K00001"]


def test_parse_provided_ko_headerless(tmp_path: Path):
    p = tmp_path / "no_header.tsv"
    p.write_text("abc K00001\ndef K00002\n")
    rows = parse_provided_ko_table(p)
    assert {r.gene_id for r in rows} == {"abc", "def"}
