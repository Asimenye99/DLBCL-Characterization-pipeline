#!/usr/bin/env python3
"""
Variant prioritisation (germline and somatic), gene-panel agnostic.

Reads VEP-annotated VCFs from two callers and retains protein-altering
variants in ANY gene. The DLBCL gene list is NOT used to filter or
prioritise: it only adds an IN_DLBCL_LIST column (YES/NO) so that genes
associated with DLBCL can be examined after annotation.

    --vcfs      tumour-only Mutect2 calls (FILTER PASS, or germline-flagged)
    --hc-vcfs   HaplotypeCaller calls after GATK hard filtering (optional)

Variants found by both callers are merged into a single row.

Each retained variant is tagged with:
    SOURCE   Mutect2, HaplotypeCaller, or Mutect2;HaplotypeCaller.
    ORIGIN   Putative_somatic or Germline_flagged (Mutect2's own inference;
             without a matched normal it is not confirmed), or Unresolved
             for variants seen only by HaplotypeCaller, which cannot
             distinguish germline from somatic.
    VAF, DP, AD   Allele fraction, depth and allele depths from the sample.
    MAX_AF   Highest population allele frequency reported by VEP.

Gene-set annotation (post-annotation, never a pathogenicity filter):
    IN_CGC / CGC_TIER / CGC_ROLE   COSMIC Cancer Gene Census membership (--census)
    IN_DLBCL_LIST                  Reddy et al. (2017) DLBCL driver gene (--flag-genes)
    GENE_SET                       Both | Census_only | Reddy_only | Neither
    Candidate variants = pathogenic-evidence variants (Level 1-3) in genes of
    the Census or the Reddy list (GENE_SET != Neither). Gene membership is not
    evidence that an individual variant is pathogenic.

Prioritisation tiers (unchanged):
    Tier 1  ClinVar pathogenic / likely pathogenic
    Tier 2  COSMIC
    Tier 3  Computational evidence (REVEL, CADD, AlphaMissense)
    Tier 4  HIGH impact consequence
    Tier 5  Other protein-altering variants
"""

import argparse
import gzip
import os
from collections import Counter


# ============================================================
# COSMIC Cancer Gene Census
# ============================================================

CENSUS = {}          # gene symbol -> {"tier": str, "role": str}
CENSUS_VIA_SYN = {}  # current symbols matched through CENSUS_ALIASES


# Verified renames between the Census symbol and the current HGNC symbol used
# by VEP. The free-text Synonyms column is NOT used: it produces false matches
# (for example FLG is listed as a synonym of FGFR1).
CENSUS_ALIASES = {"CARS1": "CARS"}


def load_census(path):
    """Read the Census CSV (columns 'Gene Symbol', 'Tier', 'Role in Cancer').
    Matching is by exact HGNC symbol plus CENSUS_ALIASES."""
    import csv
    census = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            g = (r.get("Gene Symbol") or "").strip()
            if g:
                census[g] = {"tier": (r.get("Tier") or "").strip(),
                             "role": (r.get("Role in Cancer") or "").strip()}
    via = {new: census[old] for new, old in CENSUS_ALIASES.items() if old in census}
    return census, via


def gene_set(gene, in_reddy):
    """Return (IN_CGC, CGC_TIER, CGC_ROLE, GENE_SET) for a gene symbol."""
    rec = CENSUS.get(gene) or CENSUS_VIA_SYN.get(gene)
    in_cgc = rec is not None
    if in_cgc and in_reddy:
        label = "Both"
    elif in_cgc:
        label = "Census_only"
    elif in_reddy:
        label = "Reddy_only"
    else:
        label = "Neither"
    return ("YES" if in_cgc else "NO",
            rec["tier"] if rec else "", rec["role"] if rec else "", label)


# ============================================================
# Protein-altering consequences
# ============================================================

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


# ============================================================
# Computational prediction thresholds
# ============================================================

REVEL_THRESHOLD = 0.50
CADD_THRESHOLD = 20.0
ALPHAMISSENSE_THRESHOLD = 0.564


# Thresholds can be overridden from the command line (see main()).


# ============================================================
# File and field helpers
# ============================================================

def open_vcf(path):
    """Open a plain-text or gzipped VCF."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def sample_name(path):
    """Extract the sample name from a VEP VCF filename."""
    return (
        os.path.basename(path)
        .replace(".vep.vcf.gz", "")
        .replace(".vep.vcf", "")
    )


def get_info(info, key):
    """Return the value of a key in a VCF INFO field, or an empty string."""
    for item in info.split(";"):
        if item.startswith(key + "="):
            return item.split("=", 1)[1]
    return ""


def parse_csq_header(header_lines):
    """Read the VEP CSQ column order from the VCF header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            if "Format: " not in line:
                raise RuntimeError(
                    "CSQ header found but VEP Format field is missing."
                )
            return (
                line.split("Format: ", 1)[1]
                .rstrip('">\n')
                .split("|")
            )
    raise RuntimeError("Unable to locate CSQ header.")


def to_float(value):
    """Convert an annotation value to float; None if missing or invalid."""
    if value is None:
        return None
    value = str(value).strip()
    if value in {"", ".", "-", "NA", "na", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


# ============================================================
# FILTER handling, origin and sample fields
# ============================================================

def filter_passes(filt, allowed_filters):
    """
    Accept a record only if every token in its FILTER column is allowed.

    FilterMutectCalls can write combined values (for example
    'germline;weak_evidence'), so the whole string is not compared.
    """
    return all(token in allowed_filters for token in filt.split(";"))


def variant_origin(filt):
    """
    Label the inferred origin from the FILTER column.

    Germline_flagged: FilterMutectCalls marked the call as probably germline.
    Putative_somatic: not flagged. Without a matched normal this is not
    confirmed somatic.
    """
    if "germline" in filt.split(";"):
        return "Germline_flagged"
    return "Putative_somatic"


def parse_sample(cols):
    """Return VAF, depth and allele depths from the first sample column."""
    if len(cols) < 10:
        return "", "", ""
    fields = dict(zip(cols[8].split(":"), cols[9].split(":")))
    af = fields.get("AF", "")
    ad = fields.get("AD", "")

    # HaplotypeCaller does not write AF; derive it from the allele depths.
    if not af and ad:
        try:
            counts = [int(x) for x in ad.split(",")]
            total = sum(counts)
            if len(counts) > 1 and total > 0:
                af = f"{sum(counts[1:]) / total:.3f}"
        except ValueError:
            pass

    return af, fields.get("DP", ""), ad


# ============================================================
# Evidence
# ============================================================

def clinvar_pathogenic(annotation):
    """
    True if ClinVar reports pathogenic or likely pathogenic.
    Conflicting classifications are excluded.
    """
    clin = (
        annotation.get("ClinVar_CLNSIG")
        or annotation.get("CLIN_SIG")
        or ""
    ).lower().strip()

    if not clin or "conflict" in clin:
        return False

    return "pathogenic" in clin


def clinvar_benign(annotation):
    """
    True if ClinVar reports benign or likely benign (and not also
    pathogenic). Expert-curated benign classifications override in silico
    predictors and HIGH-impact consequence, so such variants are not given
    Level 2 or Level 3 evidence. Conflicting classifications are not benign.
    """
    clin = (
        annotation.get("ClinVar_CLNSIG")
        or annotation.get("CLIN_SIG")
        or ""
    ).lower().strip()

    if not clin or "conflict" in clin or "pathogenic" in clin:
        return False

    return "benign" in clin


def evidence_level(evidence):
    """Evidence ladder: 1 = ClinVar P/LP; 2 = HIGH-impact LoF without
    ClinVar P/LP; 3 = predictor-only; '' = none."""
    if "ClinVar_P_LP" in evidence:
        return "1"
    if "HIGH_impact_LoF" in evidence:
        return "2"
    if evidence:
        return "3"
    return ""


def cosmic_variant(annotation):
    """Detect COSMIC identifiers (COSM/COSV) in Existing_variation."""
    existing = annotation.get("Existing_variation", "").upper()
    return "COSM" in existing or "COSV" in existing


def computational_evidence(annotation):
    """
    Evaluate REVEL, CADD and AlphaMissense.

    Returns (reasons, strong); strong is True if at least one criterion
    is met.
    """
    reasons = []

    revel = to_float(annotation.get("REVEL"))
    if revel is not None and revel >= REVEL_THRESHOLD:
        reasons.append(f"REVEL>={REVEL_THRESHOLD}")

    cadd = to_float(annotation.get("CADD_PHRED"))
    if cadd is not None and cadd >= CADD_THRESHOLD:
        reasons.append(f"CADD_PHRED>={CADD_THRESHOLD}")

    alpha_score = to_float(annotation.get("am_pathogenicity"))
    if alpha_score is not None and alpha_score >= ALPHAMISSENSE_THRESHOLD:
        reasons.append(f"AlphaMissense>={ALPHAMISSENSE_THRESHOLD}")

    alpha_class = annotation.get("am_class", "").strip().lower()
    if alpha_class in {"pathogenic", "likely pathogenic"}:
        reasons.append("AlphaMissense_class=" + alpha_class)

    return reasons, len(reasons) > 0


def pathogenic_evidence(annotation, min_predictors=1):
    """
    Independent evidence of pathogenicity, regardless of tier order.

    Returns a list containing any of:
        ClinVar_P_LP                       ClinVar pathogenic / likely pathogenic
        Predicted_pathogenic:<tools>       at least `min_predictors` of REVEL,
                                           CADD and AlphaMissense exceed their
                                           thresholds
        HIGH_impact_LoF                    HIGH-impact consequence (stop gain,
                                           frameshift, canonical splice site,
                                           start/stop loss)

    COSMIC membership alone is not evidence of pathogenicity and is ignored.
    """
    evidence = []

    if clinvar_benign(annotation):
        return evidence

    if clinvar_pathogenic(annotation):
        evidence.append("ClinVar_P_LP")

    reasons, _ = computational_evidence(annotation)
    predictors = set()
    for reason in reasons:
        if reason.startswith("AlphaMissense"):
            predictors.add("AlphaMissense")
        else:
            predictors.add(reason.split(">=")[0])

    if predictors and len(predictors) >= min_predictors:
        evidence.append("Predicted_pathogenic:" + "+".join(sorted(predictors)))

    if annotation.get("IMPACT", "").upper() == "HIGH":
        evidence.append("HIGH_impact_LoF")

    return evidence


def classify(annotation):
    """
    Assign a priority tier (1 highest, 5 lowest).

    Tier 5 variants are retained as candidates but do not carry strong
    pathogenic evidence. Returns (None, None) if not protein-altering.
    """
    consequence = annotation.get("Consequence", "")
    impact = annotation.get("IMPACT", "").upper()

    if clinvar_pathogenic(annotation):
        return 1, "ClinVar_pathogenic_or_likely_pathogenic"

    if cosmic_variant(annotation):
        return 2, "COSMIC"

    reasons, strong = computational_evidence(annotation)
    if strong:
        return 3, ";".join(reasons)

    if impact == "HIGH":
        return 4, "HIGH_impact"

    if any(t in PROTEIN_CONSEQUENCES for t in consequence.split("&")):
        return 5, "Protein_altering_variant_without_strong_evidence"

    return None, None


# ============================================================
# Process one VCF
# ============================================================

def process_vcf(vcf, flag_genes, allowed_filters, max_af_cutoff, max_tier,
                source="Mutect2", pathogenic_only=False, min_predictors=1):

    rows = []
    sample = sample_name(vcf)
    header_lines = []
    csq_fields = None
    cons_idx = None
    stats = Counter()

    with open_vcf(vcf) as fh:

        for line in fh:

            if line.startswith("##"):
                header_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                csq_fields = parse_csq_header(header_lines)
                cons_idx = (
                    csq_fields.index("Consequence")
                    if "Consequence" in csq_fields else None
                )
                continue

            if line.startswith("#"):
                continue

            if csq_fields is None:
                raise RuntimeError(f"CSQ header was not found in {vcf}")

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]

            stats["records"] += 1

            # ---------------------------------------------------------
            # FILTER: every token must be allowed (PASS, ., germline)
            # ---------------------------------------------------------
            if not filter_passes(filt, allowed_filters):
                stats["filtered_out"] += 1
                continue

            if source == "HaplotypeCaller":
                origin = "Unresolved"
            else:
                origin = variant_origin(filt)
            vaf, depth, allele_depth = parse_sample(cols)

            csq = get_info(info, "CSQ")
            if not csq:
                continue

            seen = set()

            for transcript in csq.split(","):

                stats["transcripts"] += 1

                values = transcript.split("|")

                # Fast pre-check: skip transcripts whose consequence is not
                # protein-altering before building the annotation dict.
                if cons_idx is not None and cons_idx < len(values):
                    if not any(
                        t in PROTEIN_CONSEQUENCES
                        for t in values[cons_idx].split("&")
                    ):
                        continue

                if len(values) < len(csq_fields):
                    values.extend([""] * (len(csq_fields) - len(values)))
                elif len(values) > len(csq_fields):
                    values = values[:len(csq_fields)]

                ann = dict(zip(csq_fields, values))

                gene = ann.get("SYMBOL", "").strip()

                if not gene:
                    # Intergenic or unannotated: nothing to model or map.
                    stats["no_gene"] += 1
                    continue

                in_list = "YES" if gene in flag_genes else "NO"
                in_cgc, cgc_tier, cgc_role, gset = gene_set(gene, in_list == "YES")
                stats["in_dlbcl_list"] += in_list == "YES"

                consequence = ann.get("Consequence", "")

                if not any(
                    t in PROTEIN_CONSEQUENCES
                    for t in consequence.split("&")
                ):
                    continue

                stats["protein"] += 1

                ensp = ann.get("ENSP", "")
                if ensp:
                    stats["ensp"] += 1

                # ---------------------------------------------------------
                # Population frequency filter. ClinVar pathogenic or
                # likely pathogenic variants are retained regardless.
                # ---------------------------------------------------------
                max_af = to_float(ann.get("MAX_AF"))

                if (
                    max_af_cutoff is not None
                    and max_af is not None
                    and max_af > max_af_cutoff
                    and not clinvar_pathogenic(ann)
                ):
                    stats["common_removed"] += 1
                    continue

                # ---------------------------------------------------------
                # Prioritisation
                # ---------------------------------------------------------
                tier, reason = classify(ann)

                if tier is None or tier > max_tier:
                    continue

                evidence = pathogenic_evidence(ann, min_predictors)

                if pathogenic_only and not evidence:
                    stats["no_pathogenic_evidence"] += 1
                    continue

                stats["classified"] += 1
                stats[f"origin_{origin}"] += 1

                # Same variant in different patients is preserved (SAMPLE).
                key = (
                    sample, chrom, pos, ref, alt, gene,
                    ann.get("HGVSp", ""), ensp,
                )

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
                    "SOURCE": source,
                    "ORIGIN": origin,
                    "VAF": vaf,
                    "DP": depth,
                    "AD": allele_depth,
                    "GENE": gene,
                    "IN_DLBCL_LIST": in_list,
                    "IN_CGC": in_cgc,
                    "CGC_TIER": cgc_tier,
                    "CGC_ROLE": cgc_role,
                    "GENE_SET": gset,
                    "EXON": ann.get("EXON", ""),
                    "DOMAINS": ann.get("DOMAINS", ""),
                    "IMPACT": ann.get("IMPACT", ""),
                    "CONSEQUENCE": consequence,
                    "HGVSc": ann.get("HGVSc", ""),
                    "HGVSp": ann.get("HGVSp", ""),
                    "Protein_position": ann.get("Protein_position", ""),
                    "Amino_acids": ann.get("Amino_acids", ""),
                    "ENSP": ensp,
                    "SWISSPROT": ann.get("SWISSPROT", ""),
                    "CANONICAL": ann.get("CANONICAL", ""),
                    "MANE": (
                        ann.get("MANE_SELECT", "")
                        or ann.get("MANE", "")
                    ),
                    "CLIN_SIG": ann.get("CLIN_SIG", ""),
                    "ClinVar_CLNSIG": ann.get("ClinVar_CLNSIG", ""),
                    "ClinVar_CLNSIGCONF": ann.get("ClinVar_CLNSIGCONF", ""),
                    "ClinVar_CLNREVSTAT": ann.get("ClinVar_CLNREVSTAT", ""),
                    "Existing_variation": ann.get("Existing_variation", ""),
                    "MAX_AF": ann.get("MAX_AF", ""),
                    "gnomADe_AF": ann.get("gnomADe_AF", ""),
                    "gnomADg_AF": ann.get("gnomADg_AF", ""),
                    "REVEL": ann.get("REVEL", ""),
                    "CADD_PHRED": ann.get("CADD_PHRED", ""),
                    "CADD_RAW": ann.get("CADD_RAW", ""),
                    "AlphaMissense_class": ann.get("am_class", ""),
                    "AlphaMissense_score": ann.get("am_pathogenicity", ""),
                    "TIER": tier,
                    "TIER_REASONS": reason,
                    "PATHOGENIC_EVIDENCE": ";".join(evidence),
                    "EVIDENCE_LEVEL": evidence_level(evidence),
                })

                stats["written"] += 1

    # ------------------------------------------------------------
    # Per-VCF summary
    # ------------------------------------------------------------
    print()
    print("====================================================")
    print(f"VCF: {vcf}  (source: {source})")
    print("====================================================")
    print("VCF records                 :", stats["records"])
    print("Removed by FILTER           :", stats["filtered_out"])
    print("VEP transcripts             :", stats["transcripts"])
    print("Transcripts without a gene  :", stats["no_gene"])
    print("In DLBCL list (flag only)   :", stats["in_dlbcl_list"])
    print("Protein-altering            :", stats["protein"])
    print("With ENSP                   :", stats["ensp"])
    print("Removed by population AF    :", stats["common_removed"])
    print("No pathogenic evidence      :", stats["no_pathogenic_evidence"])
    print()
    print("Retained by origin:")
    print("  Putative_somatic          :", stats["origin_Putative_somatic"])
    print("  Germline_flagged          :", stats["origin_Germline_flagged"])
    print("  Unresolved (HaplotypeCaller):", stats["origin_Unresolved"])
    print()
    print("Final classification:")
    for tier in range(1, 6):
        print(
            f"  Tier {tier}                    :",
            sum(1 for row in rows if row["TIER"] == tier),
        )
    print()
    print("Rows written                :", stats["written"])
    print("====================================================")
    print()

    return rows


# ============================================================
# Merge Mutect2 and HaplotypeCaller rows
# ============================================================

def row_key(row):
    return (
        row["SAMPLE"], row["CHROM"], row["POS"], row["REF"], row["ALT"],
        row["GENE"], row["HGVSp"], row["ENSP"],
    )


def merge_rows(mutect_rows, hc_rows):
    """
    One row per variant per sample.

    A variant found by both callers keeps the Mutect2 row (which carries
    the origin inference and Mutect2 allele fraction) and is labelled
    'Mutect2;HaplotypeCaller'. Variants found only by HaplotypeCaller are
    kept with ORIGIN 'Unresolved'.
    """
    merged = {}

    for row in mutect_rows:
        merged.setdefault(row_key(row), row)

    for row in hc_rows:
        key = row_key(row)
        if key in merged:
            merged[key]["SOURCE"] = "Mutect2;HaplotypeCaller"
        else:
            merged[key] = row

    return list(merged.values())


# ============================================================
# Main
# ============================================================

OUTPUT_COLUMNS = [
    "SAMPLE", "CHROM", "POS", "REF", "ALT", "FILTER",
    "SOURCE", "ORIGIN", "VAF", "DP", "AD",
    "GENE", "IMPACT", "CONSEQUENCE", "HGVSc", "HGVSp",
    "Protein_position", "Amino_acids", "ENSP", "CANONICAL", "MANE",
    "CLIN_SIG", "ClinVar_CLNSIG", "ClinVar_CLNSIGCONF",
    "ClinVar_CLNREVSTAT",
    "Existing_variation",
    "MAX_AF", "gnomADe_AF", "gnomADg_AF",
    "REVEL", "CADD_PHRED", "CADD_RAW",
    "AlphaMissense_class", "AlphaMissense_score",
    "TIER", "TIER_REASONS", "PATHOGENIC_EVIDENCE",
    # New columns are appended so that existing column positions are unchanged.
    "IN_DLBCL_LIST", "SWISSPROT", "EVIDENCE_LEVEL",
    "IN_CGC", "CGC_TIER", "CGC_ROLE", "GENE_SET", "EXON", "DOMAINS",
]


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prioritise protein-altering variants (germline-flagged and "
            "putative somatic) in any gene using ClinVar, "
            "COSMIC, CADD, REVEL and AlphaMissense."
        )
    )

    parser.add_argument(
        "--vcfs", nargs="+", required=True,
        help="VEP-annotated VCF files",
    )

    parser.add_argument(
        "--output", required=True,
        help="Output TSV file",
    )

    parser.add_argument(
        "--flag-genes", "--driver-genes", dest="flag_genes", default="",
        help=(
            "Comma-separated gene list used ONLY to fill the IN_DLBCL_LIST "
            "column. It never filters or prioritises variants. "
            "(--driver-genes is accepted as an alias.)"
        ),
    )

    parser.add_argument(
        "--pathogenic-output", default=None,
        help=(
            "Optional second TSV containing only rows with pathogenic "
            "evidence (PATHOGENIC_EVIDENCE not empty)."
        ),
    )

    parser.add_argument("--cadd-phred", type=float, default=None,
                        help="CADD PHRED threshold (default 20)")
    parser.add_argument("--revel", type=float, default=None,
                        help="REVEL threshold (default 0.5)")
    parser.add_argument("--alphamissense", type=float, default=None,
                        help="AlphaMissense threshold (default 0.564)")

    parser.add_argument(
        "--include-filters", nargs="+",
        default=["PASS", ".", "germline"],
        help=(
            "VCF FILTER tokens to accept. A record is kept only if every "
            "token in its FILTER column is in this list. Default: "
            "PASS . germline (putative somatic and germline-flagged)."
        ),
    )

    parser.add_argument(
        "--hc-vcfs", nargs="*", default=[],
        help=(
            "VEP-annotated HaplotypeCaller VCFs (hard-filtered). Optional; "
            "results are merged with the Mutect2 calls."
        ),
    )

    parser.add_argument(
        "--hc-include-filters", nargs="+", default=["PASS", "."],
        help="FILTER tokens accepted for HaplotypeCaller VCFs. Default: PASS .",
    )

    parser.add_argument(
        "--max-af", type=float, default=None,
        help=(
            "Discard variants whose VEP MAX_AF exceeds this value, unless "
            "ClinVar reports pathogenic or likely pathogenic. Variants "
            "with no frequency are kept. Default: no frequency filter."
        ),
    )

    parser.add_argument(
        "--pathogenic-only", action="store_true",
        help=(
            "Keep only variants with pathogenic evidence: ClinVar "
            "pathogenic/likely pathogenic, predicted pathogenic "
            "(see --min-predictors), or HIGH-impact loss of function. "
            "COSMIC membership alone is not sufficient."
        ),
    )

    parser.add_argument(
        "--min-predictors", type=int, default=1, choices=[1, 2, 3],
        help=(
            "Number of concordant predictors (REVEL, CADD, AlphaMissense) "
            "required for 'predicted pathogenic'. Default 1; 2 is stricter."
        ),
    )

    parser.add_argument(
        "--census", default=None,
        help="COSMIC Cancer Gene Census CSV (adds IN_CGC, CGC_TIER, CGC_ROLE, "
             "GENE_SET). Annotation only; never filters variants.",
    )
    parser.add_argument(
        "--modelling-census-tiers", nargs="+", default=["1"],
        help="Census tiers whose Level 3 missense variants go to FoldX "
             "(Level 1 missense variants in any candidate gene always do; "
             "Level 3 missense variants in Reddy genes always do). Default: 1.",
    )
    parser.add_argument(
        "--candidate-levels", nargs="+", default=["1", "2", "3"],
        help="Evidence levels that define the primary candidate set used for "
             "candidate, truncating, modelling and pathway outputs. Use '1' to "
             "restrict the analysis to ClinVar pathogenic/likely pathogenic "
             "variants. pathogenic_variants.tsv always keeps levels 1-3.",
    )
    parser.add_argument(
        "--candidate-output", default=None,
        help="TSV of candidate variants: Level 1-3 evidence in a Census or "
             "Reddy gene (GENE_SET != Neither).",
    )
    parser.add_argument(
        "--truncating-output", default=None,
        help="TSV of truncating candidates (stop_gained, frameshift, splice "
             "donor/acceptor, start/stop lost) for separate, non-FoldX "
             "characterisation.",
    )

    parser.add_argument(
        "--modelling-output", default=None,
        help=(
            "Optional TSV of missense variants (DLBCL-list genes) eligible for structural "
            "modelling: all Level 1 missense variants, plus Level 3 "
            "missense variants with >= --min-predictors concordant "
            "predictors in genes flagged IN_DLBCL_LIST (a post-annotation "
            "interpretation step; the list never filters the main tables)."
        ),
    )

    parser.add_argument(
        "--pathway-output", default=None,
        help=(
            "Optional TSV for pathway analysis: pathogenic-evidence rows in "
            "genes altered in at least --pathway-min-samples samples."
        ),
    )

    parser.add_argument(
        "--pathway-min-samples", type=int, default=1,
        help="Minimum number of samples per gene for --pathway-output.",
    )

    parser.add_argument(
        "--max-tier", type=int, default=5, choices=[1, 2, 3, 4, 5],
        help=(
            "Retain variants up to this tier. Default 5 keeps all "
            "protein-altering candidates; 4 drops Tier 5 variants that "
            "lack any pathogenic or functional evidence."
        ),
    )

    args = parser.parse_args()

    global REVEL_THRESHOLD, CADD_THRESHOLD, ALPHAMISSENSE_THRESHOLD
    if args.cadd_phred is not None:
        CADD_THRESHOLD = args.cadd_phred
    if args.revel is not None:
        REVEL_THRESHOLD = args.revel
    if args.alphamissense is not None:
        ALPHAMISSENSE_THRESHOLD = args.alphamissense

    global CENSUS, CENSUS_VIA_SYN
    if args.census:
        CENSUS, CENSUS_VIA_SYN = load_census(args.census)
        print("Census genes loaded       :", len(CENSUS),
              f"(+{len(CENSUS_VIA_SYN)} aliases)")

    allowed_filters = set(args.include_filters)
    hc_allowed_filters = set(args.hc_include_filters)

    flag_genes = {
        gene.strip()
        for gene in args.flag_genes.split(",")
        if gene.strip()
    }

    print()
    print("====================================================")
    print("VARIANT PRIORITISATION (GENE-PANEL AGNOSTIC)")
    print("====================================================")
    print()
    print("Genes in DLBCL list (flag only):", len(flag_genes))
    print()
    print("Computational thresholds:")
    print(f"REVEL >= {REVEL_THRESHOLD}")
    print(f"CADD PHRED >= {CADD_THRESHOLD}")
    print(f"AlphaMissense >= {ALPHAMISSENSE_THRESHOLD}")
    print()
    print("Accepted Mutect2 FILTER tokens:", sorted(allowed_filters))
    print("HaplotypeCaller VCFs      :", len(args.hc_vcfs))
    print("Accepted HC FILTER tokens :", sorted(hc_allowed_filters))
    print("Population MAX_AF cutoff  :", args.max_af)
    print("Maximum tier retained     :", args.max_tier)
    print("Pathogenic evidence only  :", args.pathogenic_only)
    print("Concordant predictors     :", args.min_predictors)
    print()
    print("====================================================")
    print()

    mutect_rows = []
    hc_rows = []

    for vcf in args.vcfs:
        print(f"Processing {vcf}")
        mutect_rows.extend(
            process_vcf(
                vcf, flag_genes, allowed_filters, args.max_af, args.max_tier,
                source="Mutect2",
                pathogenic_only=args.pathogenic_only,
                min_predictors=args.min_predictors,
            )
        )

    for vcf in args.hc_vcfs:
        print(f"Processing {vcf}")
        hc_rows.extend(
            process_vcf(
                vcf, flag_genes, hc_allowed_filters, args.max_af, args.max_tier,
                source="HaplotypeCaller",
                pathogenic_only=args.pathogenic_only,
                min_predictors=args.min_predictors,
            )
        )

    results = merge_rows(mutect_rows, hc_rows)

    results.sort(
        key=lambda x: (
            int(x["TIER"]),
            x["GENE"],
            x["CHROM"],
            int(x["POS"]),
        )
    )

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(args.output, "w") as out:
        out.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in results:
            out.write(
                "\t".join(str(row.get(c, "")) for c in OUTPUT_COLUMNS)
                + "\n"
            )

    if args.pathogenic_output:
        path_rows = [r for r in results if r.get("PATHOGENIC_EVIDENCE")]
        pdir = os.path.dirname(args.pathogenic_output)
        if pdir:
            os.makedirs(pdir, exist_ok=True)
        with open(args.pathogenic_output, "w") as out:
            out.write("\t".join(OUTPUT_COLUMNS) + "\n")
            for row in path_rows:
                out.write(
                    "\t".join(str(row.get(c, "")) for c in OUTPUT_COLUMNS)
                    + "\n"
                )
        def write_rows(path, rows):
            d = os.path.dirname(path)
            if d:
                os.makedirs(d, exist_ok=True)
            with open(path, "w") as fh:
                fh.write("\t".join(OUTPUT_COLUMNS) + "\n")
                for r in rows:
                    fh.write(
                        "\t".join(str(r.get(c, "")) for c in OUTPUT_COLUMNS)
                        + "\n"
                    )

        def is_missense(r):
            return "missense_variant" in r["CONSEQUENCE"].split("&")

        if args.modelling_output:
            # Post-annotation interpretation: missense variants with Level 1
            # or Level 3 evidence in DLBCL-associated genes (IN_DLBCL_LIST).
            # Level 3 rows already satisfy --min-predictors (applied in
            # pathogenic_evidence). Level 2 (loss of function) cannot be
            # modelled by FoldX.
            tiers = set(args.modelling_census_tiers)
            levels = set(args.candidate_levels)
            model_rows = [
                r for r in path_rows
                if r["EVIDENCE_LEVEL"] in levels
                and is_missense(r)
                and r["GENE_SET"] != "Neither"
                and (
                    r["EVIDENCE_LEVEL"] == "1"
                    or (r["EVIDENCE_LEVEL"] == "3"
                        and (r["IN_DLBCL_LIST"] == "YES"
                             or r["CGC_TIER"] in tiers))
                )
            ]
            write_rows(args.modelling_output, model_rows)
            print("Modelling-eligible rows   :", len(model_rows))

        if args.pathway_output:
            gene_samples = {}
            for r in path_rows:
                gene_samples.setdefault(r["GENE"], set()).add(r["SAMPLE"])
            keep = {
                g for g, ss in gene_samples.items()
                if g and len(ss) >= args.pathway_min_samples
            }
            # Post-annotation interpretation: DLBCL-associated genes only
            pw_rows = [
                r for r in path_rows
                if r["GENE"] in keep and r["GENE_SET"] != "Neither"
                and r["EVIDENCE_LEVEL"] in set(args.candidate_levels)
            ]
            write_rows(args.pathway_output, pw_rows)
            print(f"Pathway genes (candidate genes, >= {args.pathway_min_samples} "
                  f"samples):", len({r["GENE"] for r in pw_rows}))

        cand_rows = [r for r in path_rows if r["GENE_SET"] != "Neither"
                     and r["EVIDENCE_LEVEL"] in set(args.candidate_levels)]
        if args.candidate_output:
            write_rows(args.candidate_output, cand_rows)
            print("Candidate rows            :", len(cand_rows))
        if args.truncating_output:
            trunc = {"stop_gained", "frameshift_variant",
                     "splice_acceptor_variant", "splice_donor_variant",
                     "start_lost", "stop_lost"}
            write_rows(args.truncating_output,
                       [r for r in cand_rows
                        if trunc & set(r["CONSEQUENCE"].split("&"))])

        print("By evidence level (rows)  :",
              {lv: sum(1 for r in path_rows if r["EVIDENCE_LEVEL"] == lv)
               for lv in ("1", "2", "3")})
        print("Pathogenic-evidence rows  :", len(path_rows))
        print("  of which in DLBCL list  :",
              sum(1 for r in path_rows if r["IN_DLBCL_LIST"] == "YES"))

    print()
    print("====================================================")
    print("FINAL PRIORITISATION SUMMARY")
    print("====================================================")
    print("Total prioritised variants:", len(results))
    print()
    print("By origin:")
    for label in ("Putative_somatic", "Germline_flagged", "Unresolved"):
        print(
            f"  {label}: ",
            sum(1 for row in results if row["ORIGIN"] == label),
        )
    print()
    print("By source:")
    for label in (
        "Mutect2", "HaplotypeCaller", "Mutect2;HaplotypeCaller"
    ):
        print(
            f"  {label}: ",
            sum(1 for row in results if row["SOURCE"] == label),
        )
    print()
    print("By tier:")
    for tier in range(1, 6):
        print(
            f"  Tier {tier}:",
            sum(1 for row in results if row["TIER"] == tier),
        )
    print()
    print("Output file:", args.output)
    print("====================================================")
    print()


if __name__ == "__main__":
    main()
