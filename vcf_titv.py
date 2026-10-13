#!/usr/bin/env python3
"""VCF SNP Ti/Tv and per-sample diploid genotype metrics (standard library only).

Ti/Tv counts eligible REF>ALT alleles at retained VCF sites, not genotype copies.
Heterozygosity and missingness are descriptive call-set statistics, not clinical
severity classifications. No reference-genome or population inference is made.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
from pathlib import Path
from typing import Any, Iterable

BASES = frozenset("ACGT")
TRANSITIONS = frozenset({("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")})
COUNT_FIELDS = (
    "transitions", "transversions", "heterozygous", "homozygous_alt",
    "homozygous_ref", "missing", "other_ploidy",
)


def _count(value: Any, label: str) -> int:
    """Reject negative, non-integral, non-finite and ambiguous count values."""
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a non-negative integer")
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{label} must be a non-negative integer") from exc
    if not (number.is_integer() and 0 <= number <= 2**53):
        raise ValueError(f"{label} must be a non-negative integer")
    return int(number)


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def calculate_metrics(
    *, transitions: int = 0, transversions: int = 0,
    heterozygous: int = 0, homozygous_alt: int = 0,
    homozygous_ref: int = 0, missing: int = 0, other_ploidy: int = 0,
) -> dict[str, Any]:
    """Calculate descriptive ratios; undefined zero-denominator ratios are null."""
    counts = {
        key: _count(value, key)
        for key, value in locals().copy().items()
    }
    called = counts["heterozygous"] + counts["homozygous_alt"] + counts["homozygous_ref"]
    total = called + counts["missing"] + counts["other_ploidy"]
    return {
        **counts,
        "called_diploid": called,
        "genotypes_total": total,
        "ti_tv_ratio": _ratio(counts["transitions"], counts["transversions"]),
        "heterozygosity_rate": _ratio(counts["heterozygous"], called),
        "het_hom_alt_ratio": _ratio(counts["heterozygous"], counts["homozygous_alt"]),
        "missing_rate": _ratio(counts["missing"], total),
    }


def parse_vcf(lines: Iterable[str], *, pass_only: bool = False) -> dict[str, Any]:
    """Stream a plain-text VCF 4.x; count SNP alleles and per-sample GT calls.

    Missing or partially missing GT is counted as missing. Fully called
    non-diploid genotypes are counted separately, not in the diploid denominator.
    Unfiltered sites are included by default; --pass-only retains FILTER=PASS.
    """
    header = None
    samples: list[str] = []
    per_sample: dict[str, dict[str, int]] = {}
    variants = transitions = transversions = snv_alleles = 0

    for line_number, raw in enumerate(lines, start=1):
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        if line.startswith("##"):
            if header is not None:
                raise ValueError(f"Line {line_number}: metadata after column header")
            continue
        if line.startswith("#"):
            if header is not None or not line.startswith("#CHROM\t"):
                raise ValueError(f"Line {line_number}: invalid or repeated VCF header")
            header = line.split("\t")
            if header[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                raise ValueError(f"Line {line_number}: malformed VCF columns")
            if len(header) > 8 and header[8] != "FORMAT":
                raise ValueError(f"Line {line_number}: expected FORMAT column")
            samples = header[9:]
            if len(set(samples)) != len(samples) or any(not s for s in samples):
                raise ValueError(f"Line {line_number}: duplicate or blank sample names")
            per_sample = {name: {key: 0 for key in COUNT_FIELDS[2:]} for name in samples}
            continue
        if header is None:
            raise ValueError(f"Line {line_number}: missing #CHROM VCF header")

        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"Line {line_number}: expected {len(header)} tab-separated columns, got {len(fields)}")
        try:
            if int(fields[1]) < 1:
                raise ValueError
        except ValueError as exc:
            raise ValueError(f"Line {line_number}: invalid POS") from exc
        if pass_only and fields[6] != "PASS":
            continue

        alt_alleles = [] if fields[4] == "." else fields[4].split(",")
        if not fields[3] or any(not alt for alt in alt_alleles):
            raise ValueError(f"Line {line_number}: invalid REF/ALT")
        ref = fields[3].upper()
        for alt in alt_alleles:
            alternate = alt.upper()
            if len(ref) == len(alternate) == 1 and ref in BASES and alternate in BASES and ref != alternate:
                snv_alleles += 1
                if (ref, alternate) in TRANSITIONS:
                    transitions += 1
                else:
                    transversions += 1
        variants += 1

        if not samples:
            continue
        format_keys = fields[8].split(":")
        gt_index = format_keys.index("GT") if "GT" in format_keys else None
        for name, sample_field in zip(samples, fields[9:]):
            counts = per_sample[name]
            values = sample_field.split(":")
            gt = values[gt_index] if gt_index is not None and gt_index < len(values) else "."
            alleles = re.split(r"[/|]", gt)
            if any(a == "." or not a for a in alleles):
                counts["missing"] += 1
                continue
            if any(not a.isascii() or not a.isdecimal() or int(a) > len(alt_alleles) for a in alleles):
                raise ValueError(f"Line {line_number}: invalid genotype allele for sample {name}")
            if len(alleles) != 2:
                counts["other_ploidy"] += 1
            elif alleles[0] != alleles[1]:
                counts["heterozygous"] += 1
            elif alleles[0] == "0":
                counts["homozygous_ref"] += 1
            else:
                counts["homozygous_alt"] += 1

    if header is None:
        raise ValueError("Missing #CHROM VCF header")
    return {
        "variants_total": variants,
        "snv_alleles": snv_alleles,
        "transitions": transitions,
        "transversions": transversions,
        "ti_tv_ratio": _ratio(transitions, transversions),
        "samples": [
            {"sample": name, **calculate_metrics(**per_sample[name])}
            for name in samples
        ],
        "pass_only": pass_only,
    }


def analyze_file(input_path: str, *, pass_only: bool = False) -> dict[str, Any]:
    """Analyze .vcf or gzip-compressed .vcf.gz without loading all rows into RAM."""
    path = _validate_path(input_path)
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return parse_vcf(handle, pass_only=pass_only)


def _validate_path(path_str: str) -> str:
    """Reject traversal segments while allowing absolute local CLI file paths."""
    if not path_str or "\0" in path_str or ".." in path_str.replace("\\", "/").split("/"):
        raise ValueError(f"Path traversal or invalid path: {path_str}")
    return os.path.normpath(path_str)


def _csv_safe(value: Any) -> Any:
    """Avoid spreadsheet formula execution in exported user-supplied text."""
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return value


def process_batch(input_csv: str, output_csv: str) -> None:
    """Read explicit count columns, write descriptive ratios to a CSV."""
    input_csv = _validate_path(input_csv)
    output_csv = _validate_path(output_csv)
    if os.path.abspath(input_csv) == os.path.abspath(output_csv):
        raise ValueError("Input and output paths must differ")
    with open(input_csv, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = list(reader.fieldnames or [])
        if len(headers) != len(set(headers)) or not {"transitions", "transversions"}.issubset(headers):
            raise ValueError("CSV requires unique headers and transitions,transversions columns")
        output_headers = headers + [
            key for key in ("ti_tv_ratio", "heterozygosity_rate", "het_hom_alt_ratio", "missing_rate")
            if key not in headers
        ]
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Row {line_number}: too many CSV columns")
            params = {key: row[key] for key in COUNT_FIELDS if key in headers}
            result = calculate_metrics(**params)
            rows.append({
                **{key: _csv_safe(value) for key, value in row.items()},
                **{key: result[key] for key in ("ti_tv_ratio", "heterozygosity_rate", "het_hom_alt_ratio", "missing_rate")},
            })
    with open(output_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_headers)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Processed {len(rows)} records -> {output_csv}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="VCF Ti/Tv and diploid genotype statistics")
    commands = parser.add_subparsers(dest="command", required=True)

    analyze = commands.add_parser("analyze", help="Analyze a .vcf or .vcf.gz")
    analyze.add_argument("-i", "--input", required=True, help="VCF input path")
    analyze.add_argument("-o", "--output", help="Optional JSON output path")
    analyze.add_argument("--pass-only", action="store_true", help="Include FILTER=PASS records only")

    single = commands.add_parser("single", help="Compute ratios from aggregate counts")
    for key in COUNT_FIELDS:
        single.add_argument("--" + key.replace("_", "-"), type=int, default=0)

    batch = commands.add_parser("batch", help="Process counts in CSV rows")
    batch.add_argument("-i", "--input", required=True)
    batch.add_argument("-o", "--output", default="results.csv")

    args = parser.parse_args(argv)
    try:
        if args.command == "analyze":
            if args.output and os.path.abspath(_validate_path(args.input)) == os.path.abspath(_validate_path(args.output)):
                raise ValueError("Input and output paths must differ")
            result = analyze_file(args.input, pass_only=args.pass_only)
            output = json.dumps(result, indent=2, allow_nan=False)
            if args.output:
                Path(_validate_path(args.output)).write_text(output + "\n", encoding="utf-8")
            else:
                print(output)
        elif args.command == "single":
            print(json.dumps(calculate_metrics(**{key: getattr(args, key) for key in COUNT_FIELDS}), indent=2))
        else:
            process_batch(args.input, args.output)
    except (ValueError, OSError, UnicodeError, EOFError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
