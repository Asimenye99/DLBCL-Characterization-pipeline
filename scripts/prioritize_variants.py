#!/usr/bin/env python3
"""
prioritize_variants.py

Parses one or more VEP-annotated VCFs (one per sample, --vcf output with a
CSQ INFO field) and produces a single tiered candidate-variant table for
DLBCL functional characterization, tagged by sample of origin.

Tiers:
  Tier 1 - ClinVar Pathogenic / Likely_pathogenic
  Tier 2 - COSMIC-flagged AND VEP HIGH impact
  Tier 3 - CADD_PHRED >= cadd-min AND REVEL >= revel-min (missense)
  Tier 4 - AlphaMissense predicted "likely_pathogenic"
  Tier 5 - HIGH impact only (no other supporting evidence)

Only variants annotated to a gene in --driver-genes are retained.
"""

import argparse
import gzip
import os
import sys


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_csq_format(header_lines):
    """Extract the CSQ field order from the VEP VCF header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            marker = "Format: "
            idx = line.find(marker)
            if idx != -1:
                fmt_str = line[idx + len(marker):].rstrip('">\n')
                return fmt_str.split("|")
    raise ValueError("Could not find CSQ format definition in VEP VCF header.")


def get_info_field(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv.split("=", 1)[1]
    return None


def classify(csq_fields, cadd_min, revel_min):
    """Given a dict of one CSQ transcript annotation, return (tier, reasons)."""
    reasons = []
    tier = None

    impact = csq_fields.get("IMPACT", "")
    clin_sig = (csq_fields.get("ClinVar_CLNSIG") or csq_fields.get("CLIN_SIG") or "").lower()
    existing_var = (csq_fields.get("Existing_variation") or "")
    cadd = csq_fields.get("CADD_PHRED", "")
    revel = csq_fields.get("REVEL", "")
    alphamissense_class = (csq_fields.get("am_class") or csq_fields.get("AlphaMissense_class") or "")

    is_cosmic = "COSM" in existing_var or "COSV" in existing_var

    if "pathogenic" in clin_sig and "conflicting" not in clin_sig:
        tier = 1
        reasons.append(f"ClinVar={clin_sig}")
    elif is_cosmic and impact == "HIGH":
        tier = 2
        reasons.append("COSMIC-flagged")
        reasons.append("VEP_impact=HIGH")
    else:
        try:
            cadd_val = float(cadd) if cadd not in ("", ".", None) else None
        except ValueError:
            cadd_val = None
        try:
            revel_val = float(revel) if revel not in ("", ".", None) else None
        except ValueError:
            revel_val = None

        if cadd_val is not None and revel_val is not None and cadd_val >= cadd_min and revel_val >= revel_min:
            tier = 3
            reasons.append(f"CADD_PHRED={cadd_val}")
            reasons.append(f"REVEL={revel_val}")
        elif "likely_pathogenic" in alphamissense_class.lower() or "pathogenic" in alphamissense_class.lower():
            tier = 4
            reasons.append(f"AlphaMissense={alphamissense_class}")
        elif impact == "HIGH":
            tier = 5
            reasons.append("VEP_impact=HIGH (no other supporting evidence)")

    return tier, ";".join(reasons)


def sample_name_from_path(vcf_path):
    """results/annotation/WBP02.vep.vcf -> WBP02"""
    base = os.path.basename(vcf_path)
    return base.split(".vep.vcf")[0].split(".vep.vcf.gz")[0]


def process_vcf(vcf_path, driver_genes, cadd_min, revel_min):
    sample = sample_name_from_path(vcf_path)
    header_lines = []
    csq_format = None
    rows = []

    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                csq_format = parse_csq_format(header_lines)
                continue
            if csq_format is None:
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]

            if filt not in ("PASS", "."):
                continue

            csq_raw = get_info_field(info, "CSQ")
            if not csq_raw:
                continue

            for transcript_csq in csq_raw.split(","):
                values = transcript_csq.split("|")
                csq_dict = dict(zip(csq_format, values))
                gene = csq_dict.get("SYMBOL", "")

                if gene not in driver_genes:
                    continue

                tier, reasons = classify(csq_dict, cadd_min, revel_min)
                if tier is None:
                    continue

                rows.append({
                    "SAMPLE": sample,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "GENE": gene,
                    "IMPACT": csq_dict.get("IMPACT", ""),
                    "CONSEQUENCE": csq_dict.get("Consequence", ""),
                    "HGVSc": csq_dict.get("HGVSc", ""),
                    "HGVSp": csq_dict.get("HGVSp", ""),
                    "Protein_position": csq_dict.get("Protein_position", ""),
                    "Amino_acids": csq_dict.get("Amino_acids", ""),
                    "ENSP": csq_dict.get("ENSP", ""),
                    "ClinVar_CLNSIG": csq_dict.get("ClinVar_CLNSIG", ""),
                    "CADD_PHRED": csq_dict.get("CADD_PHRED", ""),
                    "REVEL": csq_dict.get("REVEL", ""),
                    "AlphaMissense_class": csq_dict.get("am_class", csq_dict.get("AlphaMissense_class", "")),
                    "Existing_variation": csq_dict.get("Existing_variation", ""),
                    "TIER": tier,
                    "TIER_REASONS": reasons,
                })

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcfs", required=True, nargs="+", help="one or more VEP-annotated VCFs (space-separated)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--driver-genes", required=True, help="comma or space separated gene symbols")
    ap.add_argument("--cadd-min", type=float, default=20)
    ap.add_argument("--revel-min", type=float, default=0.5)
    args = ap.parse_args()

    if "," in args.driver_genes:
        driver_genes = set(g.strip() for g in args.driver_genes.split(","))
    else:
        driver_genes = set(args.driver_genes.split())

    all_rows = []
    for vcf_path in args.vcfs:
        all_rows.extend(process_vcf(vcf_path, driver_genes, args.cadd_min, args.revel_min))

    all_rows.sort(key=lambda r: (r["TIER"], r["GENE"], r["SAMPLE"]))

    cols = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "GENE", "IMPACT", "CONSEQUENCE",
            "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "ENSP",
            "ClinVar_CLNSIG", "CADD_PHRED", "REVEL", "AlphaMissense_class",
            "Existing_variation", "TIER", "TIER_REASONS"]

    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in all_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    sys.stderr.write(f"Wrote {len(all_rows)} prioritized variants from {len(args.vcfs)} sample(s) to {args.output}\n")


if __name__ == "__main__":
    main()
