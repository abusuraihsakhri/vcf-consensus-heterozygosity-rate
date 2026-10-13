/* Small, dependency-free VCF 4.x statistics engine for local browser use. */
(function (root) {
  "use strict";
  const transitions = new Set(["AG", "GA", "CT", "TC"]);
  const bases = new Set(["A", "C", "G", "T"]);
  const ratio = (numerator, denominator) => denominator ? numerator / denominator : null;

  function analyze(text, options = {}) {
    if (typeof text !== "string") throw new TypeError("VCF must be text");
    const passOnly = Boolean(options.passOnly);
    const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
    let header = null;
    let names = [];
    const stats = [];
    let variants = 0, snvs = 0, ti = 0, tv = 0;

    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].replace(/\r$/, "");
      const lineNumber = i + 1;
      if (!line.trim()) continue;
      if (line.startsWith("##")) {
        if (header) throw new Error("Line " + lineNumber + ": metadata after header");
        continue;
      }
      if (line.startsWith("#")) {
        if (header || !line.startsWith("#CHROM\t")) {
          throw new Error("Line " + lineNumber + ": invalid or repeated VCF header");
        }
        header = line.split("\t");
        const required = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
        if (required.some((v, j) => header[j] !== v)) {
          throw new Error("Line " + lineNumber + ": malformed VCF columns");
        }
        if (header.length > 8 && header[8] !== "FORMAT") {
          throw new Error("Line " + lineNumber + ": expected FORMAT");
        }
        names = header.slice(9);
        if (new Set(names).size !== names.length || names.some(x => !x)) {
          throw new Error("Line " + lineNumber + ": duplicate or blank sample names");
        }
        names.forEach(name => stats.push({
          sample: name, heterozygous: 0, homozygous_alt: 0,
          homozygous_ref: 0, missing: 0, other_ploidy: 0
        }));
        continue;
      }
      if (!header) throw new Error("Line " + lineNumber + ": missing #CHROM VCF header");
      const row = line.split("\t");
      if (row.length !== header.length) {
        throw new Error("Line " + lineNumber + ": expected " + header.length + " columns, got " + row.length);
      }
      if (!/^[0-9]+$/.test(row[1]) || !Number.isSafeInteger(Number(row[1])) || Number(row[1]) < 1) {
        throw new Error("Line " + lineNumber + ": invalid POS");
      }
      if (passOnly && row[6] !== "PASS") continue;
      const ref = row[3].toUpperCase();
      const alts = row[4] === "." ? [] : row[4].split(",");
      if (!ref || alts.some(a => !a)) {
        throw new Error("Line " + lineNumber + ": invalid REF/ALT");
      }
      for (const value of alts) {
        const alt = value.toUpperCase();
        if (ref.length === 1 && alt.length === 1 && bases.has(ref) && bases.has(alt) && ref !== alt) {
          snvs++;
          if (transitions.has(ref + alt)) ti++;
          else tv++;
        }
      }
      variants++;

      if (names.length === 0) continue;
      const keys = row[8].split(":");
      const gtIndex = keys.indexOf("GT");
      for (let s = 0; s < names.length; s++) {
        const entry = stats[s];
        const values = row[s + 9].split(":");
        const gt = gtIndex >= 0 && gtIndex < values.length ? values[gtIndex] : ".";
        const alleles = gt.split(/[\/|]/);
        if (alleles.some(a => a === "." || a === "")) {
          entry.missing++;
          continue;
        }
        if (alleles.some(a => !/^[0-9]+$/.test(a) || !Number.isSafeInteger(Number(a)) || Number(a) > alts.length)) {
          throw new Error("Line " + lineNumber + ": invalid GT for sample " + names[s]);
        }
        if (alleles.length !== 2) entry.other_ploidy++;
        else if (alleles[0] !== alleles[1]) entry.heterozygous++;
        else if (Number(alleles[0]) === 0) entry.homozygous_ref++;
        else entry.homozygous_alt++;
      }
    }
    if (!header) throw new Error("Missing #CHROM VCF header");
    return {
      variants_total: variants, snv_alleles: snvs,
      transitions: ti, transversions: tv, ti_tv_ratio: ratio(ti, tv),
      samples: stats.map(s => {
        const called = s.heterozygous + s.homozygous_alt + s.homozygous_ref;
        const total = called + s.missing + s.other_ploidy;
        return Object.assign({}, s, {
          called_diploid: called, genotypes_total: total,
          heterozygosity_rate: ratio(s.heterozygous, called),
          het_hom_alt_ratio: ratio(s.heterozygous, s.homozygous_alt),
          missing_rate: ratio(s.missing, total)
        });
      }),
      pass_only: passOnly
    };
  }
  const analyzer = Object.freeze({ analyze });
  if (typeof module !== "undefined" && module.exports) module.exports = analyzer;
  root.VCFAnalyzer = analyzer;
})(typeof globalThis !== "undefined" ? globalThis : this);
