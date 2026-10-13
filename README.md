# VCF Consensus Heterozygosity Rate

A local VCF 4.x quality-control utility for variant transition/transversion (Ti/Tv) ratios and per-sample genotype statistics. It provides a dependency-free Python command-line analyzer and a separate static browser interface.

## What it measures

- **Ti/Tv:** Each eligible A/C/G/T single-base REF→ALT substitution is counted once, including each eligible ALT at multiallelic sites. Transitions are A↔G or C↔T; the remaining substitutions are transversions.
- **Per-sample calls:** Heterozygous, homozygous non-reference, homozygous reference, missing/partially missing, and fully called non-diploid genotypes (GT field).
- **Heterozygosity rate:** Heterozygous / all fully called diploid genotypes at included variant records. This is **not** population heterozygosity or a genome-wide estimate.
- **Het/Hom ALT:** Heterozygous / homozygous non-reference diploid calls.
- **Missing rate:** Missing or partially missing GT / genotype records included for that sample.

Undefined ratios (zero denominators) are **null** in JSON and **N/A** in the browser. Ti/Tv counts SNP ALT alleles, while the genotype denominators include retained SNV and non-SNV records; these metrics intentionally have different counting units.

By default all FILTER values are included; --pass-only includes exactly FILTER=PASS. Genotype quality (GQ, DP) is not filtered. Symbolic alleles and indels do not contribute to Ti/Tv. Only fully called diploid GTs contribute to the heterozygosity denominator. Results are descriptive call-set quality control, **not clinical risk assessments**.

## Browser application

The static interface is in [web/](web/). It works without a backend. Select a plain-text .vcf file (maximum 20 MB), paste VCF text, or load the example. It displays counts, per-sample statistics, and downloadable JSON/CSV. No data is transmitted by the application, stored in browser storage, or sent to an API. A GitHub Pages workflow is provided; see **Deployments** for its published status.

For local use, open web/index.html or run python -m http.server --directory web 8001 and visit http://localhost:8001.

## Python CLI

Python **3.10+** is supported. The core vcf_titv.py requires only the standard library.

~~~bash
python vcf_titv.py analyze -i sample.vcf
python vcf_titv.py analyze -i sample.vcf -o summary.json --pass-only
python vcf_titv.py analyze -i data.vcf.gz -o summary.json
python vcf_titv.py single --transitions 12 --transversions 4 --heterozygous 6 --homozygous-ref 4
python vcf_titv.py batch -i sample.csv -o results.csv
~~~

The batch CSV command expects columns transitions and transversions; optional fields are heterozygous, homozygous_alt, homozygous_ref, missing, other_ploidy and metadata. It outputs Ti/Tv, heterozygosity, Het/Hom ALT and missingness ratios, escaping formula-like spreadsheet cells. Batch CSV input contains aggregate counts, **not raw variant records**.

Install the project and tests:

~~~bash
python -m pip install -e ".[dev]"
vcf-consensus analyze -i sample.vcf
python -m pytest -q
node --test web/test_vcf.cjs
~~~

vcf-consensus is the VCF analyzer entry point. The older generic task-audit CLI remains separately available as vcf-agent or python cli.py. Its worker thresholds, mock chat and FastAPI /api/audit endpoint are **not** part of the VCF statistics calculation and must not be interpreted as validated genomic or clinical decisions.

## Assumptions and limitations

Text VCF 4.x and gzip-compressed VCF are supported by Python; the browser supports **plain-text VCF only**. Binary BCF, genotype likelihoods, phasing-aware haplotype statistics, allele-frequency estimation, reference-block handling and population-specific interpretation are outside the scope. Genotypes are counted at included records, not deduplicated across overlapping variants. Non-ACGT alleles are excluded from Ti/Tv. The Python parser streams records rather than loading the entire VCF into memory.

## Privacy, security and deployment

The browser processes files locally and has no network calls or external scripts. Downloads remain in the visitor's browser. Do not publish real participant identifiers in fixtures.

The optional agent REST server is **not** needed for browser analysis. Its audit log uses an HMAC-SHA256 chain. Set a unique, high-entropy AUDIT_SECRET_KEY environment variable to ensure consistency for a running instance. Docker Compose requires that variable:

~~~bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
docker compose up --build
~~~

Static GitHub Pages requires no credentials. The legacy agent identifier-pattern guard is not a comprehensive PHI detection or regulatory-compliance solution; restrict access to the optional API appropriately.

## Tests and source

GitHub Actions runs Python 3.10–3.13 installation, compilation, and pytest, plus Node.js syntax and parser tests. Synthetic fixtures cover multiallelic SNPs, filtered records, missing genotypes, non-diploid calls, and malformed VCF.

Core files: vcf_titv.py (Python parser and CLI), web/vcf.js (browser parser), web/app.js (UI), tests/ (regressions). The separate legacy task-audit framework remains in agents/, cli.py, enrichment.py and simulator.py.

## License

[MIT](LICENSE).
