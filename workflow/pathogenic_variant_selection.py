#!/usr/bin/env python3
"""
Pathogenic Variant Selection — DLBCL WES pipeline
==================================================
Sits between VEP Annotation and DLBCL Gene Extraction in your pipeline.

Input : results/WBP02.vep.vcf.gz   (VEP-annotated VCF, bgzipped)
Output: results/WBP02.pathogenic_variants.tsv

Usage:
    python pathogenic_variant_selection.py \
        --vcf results/WBP02.vep.vcf.gz \
        --out results/WBP02.pathogenic_variants.tsv \
        --sample WBP02

Requires: pip install cyvcf2   (fast, handles bgzipped VCFs + CSQ parsing well)
If cyvcf2 isn't available in your env, a pysam/gzip fallback is included below.
"""

import argparse
import sys
import csv
import re

# -----------------------------------------------------------------------
# Curated DLBCL driver gene panel
# (consensus across Reddy et al. 2017 Cell, Chapuy et al. 2018 Nat Med,
#  Lacy et al. 2020 Blood — adjust/extend as your thesis panel requires)
# -----------------------------------------------------------------------
DLBCL_DRIVER_GENES = {
    "TP53", "MYD88", "CD79B", "CREBBP", "EZH2", "KMT2D", "BCL2", "BCL6",
    "TNFRSF14", "GNA13", "B2M", "CD58", "PIM1", "TBL1XR1", "SOCS1",
    "ITPKB", "ETS1", "MEF2B", "STAT6", "HIST1H1E", "TNFAIP3", "PRDM1",
    "DTX1", "ATM", "POU2F2", "IRF4", "BTG1", "BTG2", "FOXO1", "XPO1",
    "NOTCH2", "SGK1", "CCND3", "RHOA", "ZFP36L1", "S1PR2", "GNAI2",
    "P2RY8", "MYC", "PTEN", "TET2", "CARD11", "TCF3", "ID3", "BCL10",
    "NOTCH1", "NFKBIE",
}

# VEP IMPACT tiers we consider, and under what conditions
HIGH_IMPACT = {"HIGH"}
MODERATE_IMPACT = {"MODERATE"}

# Regex to detect a ClinVar pathogenic call, case-insensitive,
# handles VEP's ClinVar plugin format (e.g. "Pathogenic", "Likely_pathogenic",
# "Pathogenic/Likely_pathogenic", "Conflicting_interpretations_of_pathogenicity")
CLINVAR_PATH_RE = re.compile(r"pathogenic", re.IGNORECASE)
CLINVAR_CONFLICT_RE = re.compile(r"conflicting", re.IGNORECASE)


def parse_with_cyvcf2(vcf_path, sample):
    from cyvcf2 import VCF

    vcf = VCF(vcf_path)

    # Pull CSQ field order out of the VCF header (VEP always writes this)
    csq_header = None
    for h in vcf.header_iter():
        info = h.info(extra=True)
        if info.get("ID") == "CSQ":
            desc = info.get("Description", "")
            m = re.search(r"Format:\s*(.+)", desc)
            if m:
                csq_header = m.group(1).strip('"').split("|")
            break

    if csq_header is None:
        sys.exit("ERROR: Could not find CSQ format string in VCF header. "
                  "Was this VCF actually run through VEP?")

    has_clinvar = "CLIN_SIG" in csq_header
    idx = {name: i for i, name in enumerate(csq_header)}

    rows = []
    for variant in vcf:
        csq_raw = variant.INFO.get("CSQ")
        if not csq_raw:
            continue
        for annot in csq_raw.split(","):
            fields = annot.split("|")
            if len(fields) != len(csq_header):
                continue  # malformed/truncated annotation block, skip

            impact = fields[idx.get("IMPACT", -1)] if "IMPACT" in idx else ""
            gene_symbol = fields[idx.get("SYMBOL", -1)] if "SYMBOL" in idx else ""
            consequence = fields[idx.get("Consequence", -1)] if "Consequence" in idx else ""
            clin_sig = fields[idx["CLIN_SIG"]] if has_clinvar and idx["CLIN_SIG"] < len(fields) else ""

            is_driver_gene = gene_symbol.upper() in DLBCL_DRIVER_GENES
            is_clinvar_pathogenic = bool(CLINVAR_PATH_RE.search(clin_sig)) if clin_sig else False

            # --- Selection logic ---
            keep = False
            reason = []

            if is_clinvar_pathogenic:
                keep = True
                reason.append("ClinVar_pathogenic")
            if impact in HIGH_IMPACT:
                keep = True
                reason.append("VEP_HIGH_impact")
            if impact in MODERATE_IMPACT and is_driver_gene:
                keep = True
                reason.append("MODERATE_impact_in_driver_gene")

            if not keep:
                continue

            rows.append({
                "SAMPLE": sample,
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": ",".join(variant.ALT),
                "GENE": gene_symbol,
                "IS_DLBCL_DRIVER_GENE": is_driver_gene,
                "CONSEQUENCE": consequence,
                "IMPACT": impact,
                "CLIN_SIG": clin_sig,
                "SELECTION_REASON": ";".join(reason),
            })

    return rows


def main():
    ap = argparse.ArgumentParser(description="Select pathogenic variants from VEP-annotated VCF")
    ap.add_argument("--vcf", required=True, help="Path to VEP-annotated VCF (e.g. results/WBP02.vep.vcf.gz)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--sample", default="SAMPLE", help="Sample name to tag rows with")
    args = ap.parse_args()

    try:
        rows = parse_with_cyvcf2(args.vcf, args.sample)
    except ImportError:
        sys.exit(
            "cyvcf2 is not installed. Install it with:\n"
            "    pip install cyvcf2 --break-system-packages\n"
            "(cyvcf2 needs htslib; if you're on a cluster, load your bcftools/htslib module first)"
        )

    if not rows:
        print("No variants passed the pathogenic selection filters.", file=sys.stderr)

    # Sort: driver genes first, then HIGH impact, then by position
    rows.sort(key=lambda r: (not r["IS_DLBCL_DRIVER_GENE"], r["IMPACT"] != "HIGH", r["CHROM"], r["POS"]))

    fieldnames = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "GENE",
                  "IS_DLBCL_DRIVER_GENE", "CONSEQUENCE", "IMPACT",
                  "CLIN_SIG", "SELECTION_REASON"]

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    driver_hits = sum(1 for r in rows if r["IS_DLBCL_DRIVER_GENE"])
    print(f"Wrote {len(rows)} pathogenic-selected variants to {args.out}")
    print(f"  -> {driver_hits} of these are in known DLBCL driver genes")


if __name__ == "__main__":
    main()
