"use strict";
(function () {
  const get = id => document.getElementById(id);
  const form = get("analyze-form");
  const file = get("vcf-file");
  const text = get("vcf-text");
  const status = get("status");
  const rows = get("sample-rows");
  const jsonButton = get("download-json");
  const csvButton = get("download-csv");
  let result = null;
  const format = n => n === null || n === undefined ? "N/A" : Number(n).toLocaleString("en-US", { maximumFractionDigits: 4 });
  const percentage = n => n === null || n === undefined ? "N/A" : (100 * n).toFixed(2) + "%";

  function setStatus(message, error) {
    status.textContent = message;
    status.classList.toggle("error", Boolean(error));
  }

  function clearResults() {
    result = null;
    for (const id of ["variants", "snvs", "ti", "tv", "ratio", "samples"]) {
      get("metric-" + id).textContent = "—";
    }
    rows.replaceChildren();
    const row = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 9;
    cell.textContent = "No analysis yet.";
    row.appendChild(cell);
    rows.appendChild(row);
    jsonButton.disabled = true;
    csvButton.disabled = true;
  }

  function showResults(data) {
    result = data;
    const summary = {
      variants: data.variants_total, snvs: data.snv_alleles, ti: data.transitions,
      tv: data.transversions, ratio: data.ti_tv_ratio, samples: data.samples.length
    };
    for (const [key, val] of Object.entries(summary)) get("metric-" + key).textContent = format(val);
    rows.replaceChildren();
    if (!data.samples.length) {
      const tr = document.createElement("tr"), td = document.createElement("td");
      td.colSpan = 9; td.textContent = "This VCF has no sample genotype columns.";
      tr.appendChild(td); rows.appendChild(tr);
    }
    for (const s of data.samples) {
      const tr = document.createElement("tr");
      const values = [s.sample, s.heterozygous, s.homozygous_alt, s.homozygous_ref,
        s.missing, s.other_ploidy, percentage(s.heterozygosity_rate),
        percentage(s.missing_rate), format(s.het_hom_alt_ratio)];
      for (const value of values) {
        const cell = document.createElement("td");
        cell.textContent = String(value);
        tr.appendChild(cell);
      }
      rows.appendChild(tr);
    }
    jsonButton.disabled = false; csvButton.disabled = false;
  }

  form.addEventListener("submit", async event => {
    event.preventDefault();
    clearResults();
    const selected = file.files[0];
    if (selected && selected.size > 20 * 1024 * 1024) {
      setStatus("File exceeds the 20 MB browser limit; use the Python CLI.", true);
      return;
    }
    if (selected && /\.gz$/i.test(selected.name)) {
      setStatus("Compressed VCF is supported by the Python CLI, not this browser interface.", true);
      return;
    }
    try {
      get("analyze-button").disabled = true;
      setStatus("Parsing VCF locally…");
      const source = selected ? await selected.text() : text.value;
      if (!source.trim()) throw new Error("Choose a VCF file or paste VCF text.");
      const data = window.VCFAnalyzer.analyze(source, { passOnly: get("pass-only").checked });
      showResults(data);
      setStatus("Analyzed " + data.variants_total + " variant records and " + data.samples.length + " samples locally.");
    } catch (error) {
      setStatus(error instanceof Error ? error.message : "VCF parsing failed.", true);
    } finally {
      get("analyze-button").disabled = false;
    }
  });

  get("load-example").addEventListener("click", () => {
    file.value = "";
    text.value = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n" +
      "1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n" +
      "1\t102\t.\tC\tA,T\t.\tPASS\t.\tGT\t1/2\t1/1\n" +
      "1\t103\t.\tG\tT\t.\tq10\t.\tGT\t./.\t0/1\n";
    clearResults();
    setStatus("Example loaded. Select Analyze VCF.");
  });
  form.addEventListener("reset", () => {
    window.setTimeout(() => {
      clearResults();
      setStatus("Choose a VCF file, paste text, or load the example.");
    }, 0);
  });
  function escapeCsv(value) {
    let output = String(value === null || value === undefined ? "" : value);
    if (/^\s*[=+\-@]/.test(output)) output = "'" + output;
    return '"' + output.replace(/"/g, '""') + '"';
  }
  function download(filename, content, type) {
    const url = URL.createObjectURL(new Blob([content], { type }));
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = filename;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    window.setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  jsonButton.addEventListener("click", () => {
    if (result) download("vcf-summary.json", JSON.stringify(result, null, 2) + "\n", "application/json");
  });
  csvButton.addEventListener("click", () => {
    if (!result) return;
    const keys = ["sample", "heterozygous", "homozygous_alt", "homozygous_ref",
      "missing", "other_ploidy", "called_diploid", "genotypes_total",
      "heterozygosity_rate", "het_hom_alt_ratio", "missing_rate"];
    const lines = [keys.map(escapeCsv).join(",")].concat(
      result.samples.map(s => keys.map(key => escapeCsv(s[key])).join(","))
    );
    download("vcf-samples.csv", lines.join("\r\n") + "\r\n", "text/csv;charset=utf-8");
  });
  clearResults();
})();
