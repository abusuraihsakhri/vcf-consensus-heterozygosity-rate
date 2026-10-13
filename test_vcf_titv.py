"""Smoke tests for the actual VCF command-line API."""
import json
from vcf_titv import calculate_metrics, main, process_batch


def test_vcf_titv_single(capsys):
    assert main(["single", "--transitions", "12", "--transversions", "4"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["ti_tv_ratio"] == 3.0
    assert "clinical_recommendation" not in result


def test_vcf_titv_batch(tmp_path):
    source = tmp_path / "counts.csv"
    target = tmp_path / "results.csv"
    source.write_text("sample,transitions,transversions,heterozygous,homozygous_ref\nS1,15,3,2,2\n", encoding="utf-8")
    process_batch(str(source), str(target))
    text = target.read_text(encoding="utf-8")
    assert "S1" in text
    assert "5.0" in text
    assert "0.5" in text
