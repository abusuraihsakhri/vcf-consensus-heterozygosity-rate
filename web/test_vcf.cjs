const test = require("node:test");
const assert = require("node:assert/strict");
const { analyze } = require("./vcf.js");
const header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
const rows = [
  "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0",
  "1\t20\t.\tC\tA,T\t.\tPASS\t.\tGT\t1|2\t1/1",
  "1\t30\t.\tG\tT\t.\tq10\t.\tGT\t./1\t0/1",
  "1\t40\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t1/1",
  "1\t50\t.\tA\tG\t.\tPASS\t.\tGT\t1\t0/1"
];

test("counts per-ALT SNP types and diploid genotypes", () => {
  const r = analyze(header + rows.join("\n"));
  assert.equal(r.variants_total, 5);
  assert.equal(r.snv_alleles, 5);
  assert.equal(r.transitions, 3);
  assert.equal(r.transversions, 2);
  assert.equal(r.ti_tv_ratio, 1.5);
  assert.equal(r.samples[0].heterozygosity_rate, 1);
  assert.equal(r.samples[0].missing_rate, 0.2);
  assert.equal(r.samples[0].other_ploidy, 1);
  assert.equal(r.samples[1].heterozygosity_rate, 0.4);
});

test("PASS-only filter excludes records from every metric", () => {
  const r = analyze(header + rows.join("\n"), { passOnly: true });
  assert.equal(r.variants_total, 4);
  assert.equal(r.ti_tv_ratio, 3);
  assert.equal(r.samples[0].missing, 0);
});

test("null ratios for undefined denominators and absent samples", () => {
  const r = analyze("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tG\t.\tPASS\t.\n");
  assert.equal(r.ti_tv_ratio, null);
  assert.deepEqual(r.samples, []);
});

test("rejects corrupted VCF columns, positions and genotype indices", () => {
  assert.throws(() => analyze("abc"), /missing #CHROM/);
  assert.throws(() => analyze(header + "1\t10"), /columns/);
  assert.throws(() => analyze(header + rows[0].replace("0/1", "0/2")), /invalid GT/);
  assert.throws(() => analyze(header + rows[0].replace("\t10\t", "\t0\t")), /invalid POS/);
});

test("handles BOM, CRLF and fully missing records", () => {
  const row = "1\t1\t.\tG\tA\t.\tPASS\t.\tGT\t./.\t.";
  const r = analyze("\uFEFF" + header.replace(/\n/g, "\r\n") + row + "\r\n");
  assert.equal(r.samples[0].missing, 1);
  assert.equal(r.samples[1].missing, 1);
  assert.equal(r.samples[0].heterozygosity_rate, null);
});
