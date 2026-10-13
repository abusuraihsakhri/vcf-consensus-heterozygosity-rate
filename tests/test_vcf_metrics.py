"""Synthetic VCF edge cases and biologically interpretable summary metrics."""
import gzip
import json
import pytest
from vcf_titv import analyze_file, calculate_metrics, main, parse_vcf

VCF = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	10	.	A	G	.	PASS	.	GT	0/1	0/0
1	20	.	C	A,T	.	PASS	.	GT	1|2	1/1
1	30	.	G	T	.	q10	.	GT	./1	0/1
1	40	.	A	AT	.	PASS	.	GT	0/1	1/1
1	50	.	A	G	.	PASS	.	GT	1	0/1
"""


def test_multiallelic_and_sample_metrics():
    result = parse_vcf(VCF.splitlines())
    assert (result["variants_total"], result["snv_alleles"]) == (5, 5)
    assert (result["transitions"], result["transversions"]) == (3, 2)
    assert result["ti_tv_ratio"] == 1.5
    first, second = result["samples"]
    assert first["heterozygous"] == 3
    assert first["other_ploidy"] == 1
    assert first["missing"] == 1
    assert first["heterozygosity_rate"] == 1
    assert first["missing_rate"] == 0.2
    assert second["homozygous_alt"] == 2
    assert second["homozygous_ref"] == 1
    assert second["heterozygosity_rate"] == 0.4


def test_pass_filter_changes_sites_and_genotypes():
    result = parse_vcf(VCF.splitlines(), pass_only=True)
    assert result["variants_total"] == 4
    assert result["ti_tv_ratio"] == 3
    assert result["samples"][0]["missing"] == 0


def test_undefined_denominators_are_null():
    result = parse_vcf(["#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
                        "1	1	.	A	G	.	PASS	."])
    assert result["ti_tv_ratio"] is None
    assert result["samples"] == []
    assert calculate_metrics()["missing_rate"] is None


@pytest.mark.parametrize("vcf", [
    "1\t1\t.\tA\tG\t.\tPASS\t.",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t0\t.\tA\tG\t.\tPASS\t.",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t2\t.\tA\tG\t.\tPASS\t.\tGT\t0/2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t.\tA\tG",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1",
])
def test_malformed_vcf_rejected(vcf):
    with pytest.raises(ValueError):
        parse_vcf(vcf.splitlines())


def test_gzip_supported(tmp_path):
    path = tmp_path / "test.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(VCF)
    assert analyze_file(str(path))["variants_total"] == 5


def test_cli_json_output(tmp_path):
    source = tmp_path / "test.vcf"
    output = tmp_path / "summary.json"
    source.write_text(VCF, encoding="utf-8")
    assert main(["analyze", "-i", str(source), "-o", str(output), "--pass-only"]) == 0
    assert json.loads(output.read_text())["variants_total"] == 4


def test_batch_formula_escaping(tmp_path):
    from vcf_titv import process_batch
    source = tmp_path / "counts.csv"
    output = tmp_path / "summary.csv"
    source.write_text("sample,transitions,transversions\n=cmd,2,1\n", encoding="utf-8")
    process_batch(str(source), str(output))
    assert "'=cmd" in output.read_text()
