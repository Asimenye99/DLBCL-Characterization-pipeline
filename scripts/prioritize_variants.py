#!/usr/bin/env python3

"""
prioritize_variants.py

Prioritize protein-altering variants from VEP-annotated VCF files
for downstream DLBCL structural modelling.

Priority:
1. ClinVar pathogenic / likely pathogenic
2. COSMIC variants
3. HIGH impact variants
4. Missense variants in DLBCL driver genes
"""

import argparse
import gzip
import os
import sys


PROTEIN_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "start_lost",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def sample_name(path):
    return (
        os.path.basename(path)
        .replace(".vep.vcf.gz", "")
        .replace(".vep.vcf", "")
    )


def get_info(info, key):
    """Extract INFO field value."""
    for item in info.split(";"):
        if item.startswith(key + "="):
            return item.split("=", 1)[1]
    return ""


def parse_csq_header(header_lines):
    """Obtain the order of VEP CSQ columns dynamically."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            fmt = (
                line.split("Format: ")[1]
                .rstrip('">\n')
                .split("|")
            )
            return fmt

    raise RuntimeError("Unable to locate CSQ header.")


def classify(annotation):
    """Assign variant priority."""
    consequence = annotation.get("Consequence", "")
    impact = annotation.get("IMPACT", "")

    clin = (
        annotation.get("ClinVar_CLNSIG")
        or annotation.get("CLIN_SIG")
        or ""
    ).lower()

    existing = annotation.get("Existing_variation", "")

    if "pathogenic" in clin and "conflicting" not in clin:
        return 1, "ClinVar_pathogenic"

    if "COSM" in existing or "COSV" in existing:
        return 2, "COSMIC"

    if impact == "HIGH":
        return 3, "HIGH_impact"

    if "missense_variant" in consequence:
        return 4, "Missense"

    return None, None


def process_vcf(vcf, drivers, allowed_filters):
    rows = []
    sample = sample_name(vcf)
    header_lines = []
    csq_fields = None

    stats = {
        "transcripts": 0,
        "driver": 0,
        "protein": 0,
        "ensp": 0,
        "classified": 0,
        "written": 0,
    }

    with open_vcf(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                csq_fields = parse_csq_header(header_lines)
                continue

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]

            if filt not in allowed_filters:
                continue

            csq = get_info(info, "CSQ")
            if not csq:
                continue

            seen = set()

            for transcript in csq.split(","):
                stats["transcripts"] += 1

                values = transcript.split("|")

                # pad missing fields
                if len(values) < len(csq_fields):
                    values.extend([""] * (len(csq_fields) - len(values)))

                ann = dict(zip(csq_fields, values))

                gene = ann.get("SYMBOL", "").strip()
                if gene not in drivers:
                    continue
                stats["driver"] += 1

                consequence = ann.get("Consequence", "")
                parts = consequence.split("&")
                matches = [t for t in parts if t in PROTEIN_CONSEQUENCES]

                protein = bool(matches)
                if not protein:
                    continue
                stats["protein"] += 1

                ensp = ann.get("ENSP", "")
                if ensp:
                    stats["ensp"] += 1

                tier, reason = classify(ann)
                if tier is None:
                    continue
                stats["classified"] += 1

                key = (chrom, pos, ref, alt, gene, ann.get("HGVSp", ""), ensp)
                if key in seen:
                    continue
                seen.add(key)

                rows.append({
                    "SAMPLE": sample,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "FILTER": filt,
                    "GENE": gene,
                    "IMPACT": ann.get("IMPACT", ""),
                    "CONSEQUENCE": consequence,
                    "HGVSc": ann.get("HGVSc", ""),
                    "HGVSp": ann.get("HGVSp", ""),
                    "Protein_position": ann.get("Protein_position", ""),
                    "Amino_acids": ann.get("Amino_acids", ""),
                    "ENSP": ensp,
                    "CANONICAL": ann.get("CANONICAL", ""),
                    "MANE": ann.get("MANE_SELECT", "") or ann.get("MANE", ""),
                    "CLIN_SIG": ann.get("CLIN_SIG", ""),
                    "ClinVar_CLNSIG": ann.get("ClinVar_CLNSIG", ""),
                    "Existing_variation": ann.get("Existing_variation", ""),
                    "TIER": tier,
                    "TIER_REASONS": reason,
                })

                stats["written"] += 1

    print("\n========== DEBUG ==========")
    print("Total transcripts :", stats["transcripts"])
    print("Driver genes      :", stats["driver"])
    print("Protein altering  :", stats["protein"])
    print("With ENSP         :", stats["ensp"])
    print("Classified        :", stats["classified"])
    print("Rows written      :", stats["written"])
    print("===========================\n")

    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcfs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--driver-genes", required=True)
    p.add_argument(
        "--include-filters",
        nargs="+",
        default=["PASS", "."],
        help=(
            "FILTER values to accept (default: PASS only). "
            "e.g. --include-filters PASS . germline  to also review "
            "germline-flagged variants."
        ),
    )
    args = p.parse_args()
    allowed_filters = set(args.include_filters)

    drivers = {
        x.strip()
        for x in args.driver_genes.replace(",", " ").split()
    }

    print("Driver genes loaded:")
    print(sorted(drivers))
    print()

    results = []
    for vcf in args.vcfs:
        print(f"Processing {vcf}")
        results.extend(process_vcf(vcf, drivers, allowed_filters))

    results.sort(
        key=lambda x: (x["TIER"], x["GENE"], x["CHROM"], int(x["POS"]))
    )

    columns = [
        "SAMPLE",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "GENE",
        "IMPACT",
        "CONSEQUENCE",
        "HGVSc",
        "HGVSp",
        "Protein_position",
        "Amino_acids",
        "ENSP",
        "ClinVar_CLNSIG",
        "Existing_variation",
        "TIER",
        "TIER_REASONS",
    ]

    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in results:
            out.write(
                "\t".join(str(row.get(c, "")) for c in columns) + "\n"
            )

    print()
    print("===================================")
    print(f"Prioritized variants : {len(results)}")
    print(f"Output file          : {args.output}")
    print("===================================")


if __name__ == "__main__":
    main()
