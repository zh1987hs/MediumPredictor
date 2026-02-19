from genome_nutrition_predictor.rules import load_rules
from genome_nutrition_predictor.scoring import index_ko_to_genes, score_pathway
from genome_nutrition_predictor.models import AnnotationRecord


def test_rules_load_contains_biotin():
    rules = load_rules()
    ids = {p["id"] for p in rules["pathways"]}
    assert "b7_de_novo" in ids


def test_pathway_scoring_present_for_biotin():
    rules = load_rules()
    biotin = [p for p in rules["pathways"] if p["id"] == "b7_de_novo"][0]
    ann = [
        AnnotationRecord(gene_id="a", kos=["K00652"]),
        AnnotationRecord(gene_id="b", kos=["K00833"]),
        AnnotationRecord(gene_id="c", kos=["K01935"]),
        AnnotationRecord(gene_id="d", kos=["K01012"]),
    ]
    idx = index_ko_to_genes(ann)
    res = score_pathway(biotin, idx)
    assert res.status == "Present"
    assert res.score_0_100 >= 95
