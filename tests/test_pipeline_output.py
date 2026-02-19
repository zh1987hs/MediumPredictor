import json
from pathlib import Path

from genome_nutrition_predictor.pipeline import run_prediction


def test_toy_output(tmp_path: Path):
    out = tmp_path / "out"
    run_prediction(
        out_dir=out,
        annotation_mode="provided_ko",
        annotations="examples/toy_ko.tsv",
        proteins=None,
        genome=None,
        rules_path=None,
        db_dir=None,
        offline=True,
        min_completeness_threshold=0.8,
        threads=1,
    )
    for name in ["summary.json", "report.md", "pathways.tsv", "evidence.tsv", "medium_suggestion.tsv", "metadata.json"]:
        assert (out / name).exists()
    summary = json.loads((out / "summary.json").read_text())
    assert "auxotrophy" in summary and "capabilities" in summary
