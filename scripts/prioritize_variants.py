#!/usr/bin/env python3

import argparse
import gzip
import os


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


# ============================================================
# Open VCF
# ============================================================

def open_vcf(path):
    """
    Open plain-text or gzipped VCF.
    """

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r")


# ============================================================
# Sample name
# ============================================================

def sample_name(path):
    """
    Extract sample name from VEP VCF filename.
    """

    return (
        os.path.basename(path)
        .replace(".vep.vcf.gz", "")
        .replace(".vep.vcf", "")
    )


# ============================================================
# Extract INFO field
# ============================================================

def get_info(info, key):
    """
    Extract a value from a VCF INFO field.

    Example:
        CSQ=...
        returns the CSQ value.
    """

    for item in info.split(";"):

        if item.startswith(key + "="):

            return item.split("=", 1)[1]

    return ""


# ============================================================
# Parse VEP CSQ header
# ============================================================

def parse_csq_header(header_lines):
    """
    Obtain the order of VEP CSQ columns dynamically
    from the VCF header.

    This prevents the script from assuming a fixed
    VEP annotation-column order.
    """

    for line in header_lines:

        if line.startswith("##INFO=<ID=CSQ"):

            if "Format: " not in line:

                raise RuntimeError(
                    "CSQ header found but VEP Format field is missing."
                )

            fmt = (
                line.split("Format: ", 1)[1]
                .rstrip('">\n')
                .split("|")
            )

            return fmt

    raise RuntimeError(
        "Unable to locate CSQ header."
    )


# ============================================================
# Safe conversion to float
# ============================================================

def to_float(value):
    """
    Convert annotation value to float.

    Returns None for missing or invalid values.
    """

    if value is None:

        return None

    value = str(value).strip()

    if value in {
        "",
        ".",
        "-",
        "NA",
        "na",
        "nan",
        "NaN",
    }:

        return None

    try:

        return float(value)

    except (ValueError, TypeError):

        return None


# ============================================================
# ClinVar pathogenic classification
# ============================================================

def clinvar_pathogenic(annotation):
    """
    Determine whether ClinVar indicates pathogenic
    or likely pathogenic classification.

    Conflicting classifications are excluded.

    Accepted examples include:

        Pathogenic
        Likely pathogenic
        likely_pathogenic
        pathogenic/likely_pathogenic
    """

    clin = (
        annotation.get("ClinVar_CLNSIG")
        or annotation.get("CLIN_SIG")
        or ""
    ).lower().strip()

    if not clin:

        return False

    # Exclude explicit conflicts
    if "conflict" in clin:

        return False

    # Explicit pathogenic terminology
    if "pathogenic" in clin:

        return True

    return False


# ============================================================
# COSMIC
# ============================================================

def cosmic_variant(annotation):
    """
    Detect COSMIC identifiers from VEP Existing_variation.

    Examples:
        COSM12345
        COSV12345
    """

    existing = annotation.get(
        "Existing_variation",
        ""
    )

    return (
        "COSM" in existing.upper()
        or "COSV" in existing.upper()
    )


# ============================================================
# Computational evidence
# ============================================================

def computational_evidence(annotation):
    """
    Evaluate REVEL, CADD and AlphaMissense.

    A variant receives computational evidence if at least
    one of the following criteria is met:

        REVEL >= 0.50
        CADD PHRED >= 20
        AlphaMissense score >= 0.564
        AlphaMissense class = pathogenic
        AlphaMissense class = likely pathogenic

    Returns:

        reasons
            List describing all computational evidence found.

        strong
            True if at least one computational criterion
            is satisfied.
    """

    reasons = []

    # ========================================================
    # REVEL
    # ========================================================

    revel = to_float(
        annotation.get("REVEL")
    )

    if (
        revel is not None
        and revel >= REVEL_THRESHOLD
    ):

        reasons.append(
            f"REVEL>={REVEL_THRESHOLD}"
        )

    # ========================================================
    # CADD PHRED
    # ========================================================

    cadd = to_float(
        annotation.get("CADD_PHRED")
    )

    if (
        cadd is not None
        and cadd >= CADD_THRESHOLD
    ):

        reasons.append(
            f"CADD_PHRED>={CADD_THRESHOLD}"
        )

    # ========================================================
    # AlphaMissense score
    # ========================================================

    alpha_score = to_float(
        annotation.get("am_pathogenicity")
    )

    if (
        alpha_score is not None
        and alpha_score >= ALPHAMISSENSE_THRESHOLD
    ):

        reasons.append(
            f"AlphaMissense>={ALPHAMISSENSE_THRESHOLD}"
        )

    # ========================================================
    # AlphaMissense class
    # ========================================================

    alpha_class = (
        annotation.get(
            "am_class",
            ""
        )
        .strip()
        .lower()
    )

    if alpha_class in {
        "pathogenic",
        "likely pathogenic",
    }:

        reasons.append(
            "AlphaMissense_class=" + alpha_class
        )

    # ========================================================
    # Final computational evidence decision
    # ========================================================

    strong = len(reasons) > 0

    return reasons, strong


# ============================================================
# Variant classification
# ============================================================

def classify(annotation):
    """
    Assign variant priority.

    Tier 1:
        ClinVar pathogenic / likely pathogenic

    Tier 2:
        COSMIC

    Tier 3:
        Strong computational evidence:
        REVEL
        CADD
        AlphaMissense

    Tier 4:
        HIGH impact consequence

    Tier 5:
        Other protein-altering variants

    IMPORTANT:
        Tier 5 variants are retained as candidates but are
        not considered strong pathogenic evidence.
    """

    consequence = annotation.get(
        "Consequence",
        ""
    )

    impact = annotation.get(
        "IMPACT",
        ""
    ).upper()

    # ========================================================
    # TIER 1 — CLINVAR
    # ========================================================

    if clinvar_pathogenic(annotation):

        return (
            1,
            "ClinVar_pathogenic_or_likely_pathogenic"
        )

    # ========================================================
    # TIER 2 — COSMIC
    # ========================================================

    if cosmic_variant(annotation):

        return (
            2,
            "COSMIC"
        )

    # ========================================================
    # TIER 3 — COMPUTATIONAL EVIDENCE
    # ========================================================

    computational_reasons, strong = (
        computational_evidence(annotation)
    )

    if strong:

        return (
            3,
            ";".join(computational_reasons)
        )

    # ========================================================
    # TIER 4 — HIGH IMPACT
    # ========================================================

    if impact == "HIGH":

        return (
            4,
            "HIGH_impact"
        )

    # ========================================================
    # TIER 5 — OTHER PROTEIN-ALTERING VARIANTS
    # ========================================================

    consequence_terms = consequence.split("&")

    if any(
        term in PROTEIN_CONSEQUENCES
        for term in consequence_terms
    ):

        return (
            5,
            "Protein_altering_variant_without_strong_evidence"
        )

    return None, None


# ============================================================
# Process VCF
# ============================================================

def process_vcf(
    vcf,
    drivers,
    allowed_filters
):

    rows = []

    sample = sample_name(vcf)

    header_lines = []
    csq_fields = None

    stats = {
        "records": 0,
        "transcripts": 0,
        "driver": 0,
        "protein": 0,
        "ensp": 0,
        "clinvar": 0,
        "cosmic": 0,
        "computational": 0,
        "high_impact": 0,
        "tier5": 0,
        "classified": 0,
        "written": 0,
    }

    with open_vcf(vcf) as fh:

        for line in fh:

            # =================================================
            # VCF metadata
            # =================================================

            if line.startswith("##"):

                header_lines.append(line)

                continue

            # =================================================
            # VCF column header
            # =================================================

            if line.startswith("#CHROM"):

                csq_fields = parse_csq_header(
                    header_lines
                )

                continue

            if line.startswith("#"):

                continue

            # =================================================
            # Make sure CSQ header was found
            # =================================================

            if csq_fields is None:

                raise RuntimeError(
                    f"CSQ header was not found in {vcf}"
                )

            # =================================================
            # VCF fields
            # =================================================

            cols = line.rstrip().split("\t")

            if len(cols) < 8:

                continue

            (
                chrom,
                pos,
                vid,
                ref,
                alt,
                qual,
                filt,
                info
            ) = cols[:8]

            stats["records"] += 1

            # =================================================
            # FILTER
            # =================================================

            if filt not in allowed_filters:

                continue

            # =================================================
            # VEP CSQ
            # =================================================

            csq = get_info(
                info,
                "CSQ"
            )

            if not csq:

                continue

            # =================================================
            # Process each transcript
            # =================================================

            seen = set()

            for transcript in csq.split(","):

                stats["transcripts"] += 1

                values = transcript.split("|")

                # ---------------------------------------------
                # Pad missing fields
                # ---------------------------------------------

                if len(values) < len(csq_fields):

                    values.extend(
                        [""] *
                        (
                            len(csq_fields)
                            - len(values)
                        )
                    )

                # ---------------------------------------------
                # If there are unexpectedly more values,
                # truncate to the CSQ header length.
                # ---------------------------------------------

                if len(values) > len(csq_fields):

                    values = values[
                        :len(csq_fields)
                    ]

                # ---------------------------------------------
                # Build annotation dictionary
                # ---------------------------------------------

                ann = dict(
                    zip(
                        csq_fields,
                        values
                    )
                )

                # =================================================
                # Gene
                # =================================================

                gene = ann.get(
                    "SYMBOL",
                    ""
                ).strip()

                if gene not in drivers:

                    continue

                stats["driver"] += 1

                # =================================================
                # Consequence
                # =================================================

                consequence = ann.get(
                    "Consequence",
                    ""
                )

                consequence_terms = (
                    consequence.split("&")
                )

                matches = [
                    term
                    for term in consequence_terms
                    if term in PROTEIN_CONSEQUENCES
                ]

                if not matches:

                    continue

                stats["protein"] += 1

                # =================================================
                # ENSP
                # =================================================

                ensp = ann.get(
                    "ENSP",
                    ""
                )

                if ensp:

                    stats["ensp"] += 1

                # =================================================
                # Extract annotation values
                # =================================================

                revel = ann.get(
                    "REVEL",
                    ""
                )

                cadd_phred = ann.get(
                    "CADD_PHRED",
                    ""
                )

                cadd_raw = ann.get(
                    "CADD_RAW",
                    ""
                )

                alpha_class = ann.get(
                    "am_class",
                    ""
                )

                alpha_score = ann.get(
                    "am_pathogenicity",
                    ""
                )

                clin_sig = ann.get(
                    "ClinVar_CLNSIG",
                    ""
                )

                existing_variation = ann.get(
                    "Existing_variation",
                    ""
                )

                # =================================================
                # Evidence counts
                # =================================================

                if clinvar_pathogenic(ann):

                    stats["clinvar"] += 1

                if cosmic_variant(ann):

                    stats["cosmic"] += 1

                computational_reasons, computational_strong = (
                    computational_evidence(ann)
                )

                if computational_strong:

                    stats["computational"] += 1

                if ann.get(
                    "IMPACT",
                    ""
                ).upper() == "HIGH":

                    stats["high_impact"] += 1

                # =================================================
                # Prioritisation
                # =================================================

                tier, reason = classify(
                    ann
                )

                if tier is None:

                    continue

                stats["classified"] += 1

                if tier == 5:

                    stats["tier5"] += 1

                # =================================================
                # Deduplication
                #
                # Include SAMPLE so the same biological variant
                # occurring in different patients is preserved.
                # =================================================

                key = (
                    sample,
                    chrom,
                    pos,
                    ref,
                    alt,
                    gene,
                    ann.get(
                        "HGVSp",
                        ""
                    ),
                    ensp,
                )

                if key in seen:

                    continue

                seen.add(key)

                # =================================================
                # Save row
                # =================================================

                rows.append({

                    "SAMPLE": sample,

                    "CHROM": chrom,

                    "POS": pos,

                    "REF": ref,

                    "ALT": alt,

                    "FILTER": filt,

                    "GENE": gene,

                    "IMPACT": ann.get(
                        "IMPACT",
                        ""
                    ),

                    "CONSEQUENCE": consequence,

                    "HGVSc": ann.get(
                        "HGVSc",
                        ""
                    ),

                    "HGVSp": ann.get(
                        "HGVSp",
                        ""
                    ),

                    "Protein_position": ann.get(
                        "Protein_position",
                        ""
                    ),

                    "Amino_acids": ann.get(
                        "Amino_acids",
                        ""
                    ),

                    "ENSP": ensp,

                    "CANONICAL": ann.get(
                        "CANONICAL",
                        ""
                    ),

                    "MANE": (
                        ann.get(
                            "MANE_SELECT",
                            ""
                        )
                        or ann.get(
                            "MANE",
                            ""
                        )
                    ),

                    # ---------------------------------------------
                    # ClinVar
                    # ---------------------------------------------

                    "CLIN_SIG": ann.get(
                        "CLIN_SIG",
                        ""
                    ),

                    "ClinVar_CLNSIG": clin_sig,

                    "ClinVar_CLNSIGCONF": ann.get(
                        "ClinVar_CLNSIGCONF",
                        ""
                    ),

                    "ClinVar_CLNREVSTAT": ann.get(
                        "ClinVar_CLNREVSTAT",
                        ""
                    ),

                    # ---------------------------------------------
                    # COSMIC / Existing variation
                    # ---------------------------------------------

                    "Existing_variation": existing_variation,

                    # ---------------------------------------------
                    # Computational predictions
                    # ---------------------------------------------

                    "REVEL": revel,

                    "CADD_PHRED": cadd_phred,

                    "CADD_RAW": cadd_raw,

                    "AlphaMissense_class": alpha_class,

                    "AlphaMissense_score": alpha_score,

                    # ---------------------------------------------
                    # Priority
                    # ---------------------------------------------

                    "TIER": tier,

                    "TIER_REASONS": reason,

                })

                stats["written"] += 1

    # ========================================================
    # Debug / statistics
    # ========================================================

    print()
    print("====================================================")
    print(f"VCF: {vcf}")
    print("====================================================")

    print(
        "VCF records              :",
        stats["records"]
    )

    print(
        "VEP transcripts          :",
        stats["transcripts"]
    )

    print(
        "Driver-gene transcripts  :",
        stats["driver"]
    )

    print(
        "Protein-altering         :",
        stats["protein"]
    )

    print(
        "With ENSP                :",
        stats["ensp"]
    )

    print()
    print("Evidence detected:")
    print(
        "ClinVar pathogenic       :",
        stats["clinvar"]
    )

    print(
        "COSMIC                   :",
        stats["cosmic"]
    )

    print(
        "Computational evidence  :",
        stats["computational"]
    )

    print(
        "HIGH impact              :",
        stats["high_impact"]
    )

    print()
    print("Final classification:")
    print(
        "Tier 1                   :",
        sum(
            1
            for row in rows
            if row["TIER"] == 1
        )
    )

    print(
        "Tier 2                   :",
        sum(
            1
            for row in rows
            if row["TIER"] == 2
        )
    )

    print(
        "Tier 3                   :",
        sum(
            1
            for row in rows
            if row["TIER"] == 3
        )
    )

    print(
        "Tier 4                   :",
        sum(
            1
            for row in rows
            if row["TIER"] == 4
        )
    )

    print(
        "Tier 5                   :",
        sum(
            1
            for row in rows
            if row["TIER"] == 5
        )
    )

    print()
    print(
        "Classified               :",
        stats["classified"]
    )

    print(
        "Rows written             :",
        stats["written"]
    )

    print("====================================================")
    print()

    return rows


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prioritise protein-altering variants in "
            "DLBCL driver genes using ClinVar, COSMIC, "
            "CADD, REVEL and AlphaMissense."
        )
    )

    # ========================================================
    # Input VCFs
    # ========================================================

    parser.add_argument(
        "--vcfs",
        nargs="+",
        required=True,
        help="VEP-annotated VCF files"
    )

    # ========================================================
    # Output
    # ========================================================

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file"
    )

    # ========================================================
    # Driver genes
    #
    # IMPORTANT:
    # Snakefile passes these as a comma-separated list.
    # ========================================================

    parser.add_argument(
        "--driver-genes",
        required=True,
        help=(
            "Comma-separated list of DLBCL driver genes"
        )
    )

    # ========================================================
    # Accepted FILTER values
    # ========================================================

    parser.add_argument(
        "--include-filters",
        nargs="+",
        default=[
            "PASS",
            "."
        ],
        help=(
            "VCF FILTER values to accept. "
            "For the current pipeline use: "
            "PASS .  and Variants filtered as germline are excluded."
        ),
    )

    args = parser.parse_args()

    # ========================================================
    # Allowed filters
    # ========================================================

    allowed_filters = set(
        args.include_filters
    )

    # ========================================================
    # Load driver genes
    #
    # The Snakefile passes:
    #
    # KMT2D,TP53,BCL6,...
    # ========================================================

    drivers = {
        gene.strip()
        for gene in args.driver_genes.split(",")
        if gene.strip()
    }

    # ========================================================
    # Print driver genes
    # ========================================================

    print()
    print("====================================================")
    print("DLBCL VARIANT PRIORITISATION")
    print("====================================================")

    print()
    print("Driver genes loaded:")

    print(
        sorted(drivers)
    )

    print(
        "Number of driver genes:",
        len(drivers)
    )

    print()
    print("Computational thresholds:")

    print(
        f"REVEL >= {REVEL_THRESHOLD}"
    )

    print(
        f"CADD PHRED >= {CADD_THRESHOLD}"
    )

    print(
        f"AlphaMissense >= {ALPHAMISSENSE_THRESHOLD}"
    )

    print()
    print(
        "Accepted VCF FILTER values:",
        sorted(allowed_filters)
    )

    print()
    print("====================================================")
    print()

    # ========================================================
    # Process VCFs
    # ========================================================

    results = []

    for vcf in args.vcfs:

        print(
            f"Processing {vcf}"
        )

        vcf_results = process_vcf(
            vcf,
            drivers,
            allowed_filters
        )

        results.extend(
            vcf_results
        )

    # ========================================================
    # Sort
    #
    # Priority:
    # Tier 1 -> Tier 5
    # then gene
    # then chromosome
    # then position
    # ========================================================

    results.sort(
        key=lambda x: (
            int(x["TIER"]),
            x["GENE"],
            x["CHROM"],
            int(x["POS"])
        )
    )

    # ========================================================
    # Output columns
    # ========================================================

    columns = [

        # ----------------------------------------------------
        # Variant identification
        # ----------------------------------------------------

        "SAMPLE",

        "CHROM",

        "POS",

        "REF",

        "ALT",

        "FILTER",

        # ----------------------------------------------------
        # Gene / consequence
        # ----------------------------------------------------

        "GENE",

        "IMPACT",

        "CONSEQUENCE",

        "HGVSc",

        "HGVSp",

        "Protein_position",

        "Amino_acids",

        "ENSP",

        "CANONICAL",

        "MANE",

        # ----------------------------------------------------
        # ClinVar
        # ----------------------------------------------------

        "CLIN_SIG",

        "ClinVar_CLNSIG",

        "ClinVar_CLNSIGCONF",

        "ClinVar_CLNREVSTAT",

        # ----------------------------------------------------
        # COSMIC / existing variation
        # ----------------------------------------------------

        "Existing_variation",

        # ----------------------------------------------------
        # Computational predictors
        # ----------------------------------------------------

        "REVEL",

        "CADD_PHRED",

        "CADD_RAW",

        "AlphaMissense_class",

        "AlphaMissense_score",

        # ----------------------------------------------------
        # Prioritisation
        # ----------------------------------------------------

        "TIER",

        "TIER_REASONS",
    ]

    # ========================================================
    # Create output directory
    # ========================================================

    output_dir = os.path.dirname(
        args.output
    )

    if output_dir:

        os.makedirs(
            output_dir,
            exist_ok=True
        )

    # ========================================================
    # Write output
    # ========================================================

    with open(
        args.output,
        "w"
    ) as out:

        # Header
        out.write(
            "\t".join(columns)
            + "\n"
        )

        # Data
        for row in results:

            out.write(
                "\t".join(
                    str(
                        row.get(
                            column,
                            ""
                        )
                    )
                    for column in columns
                )
                + "\n"
            )

    # ========================================================
    # Final summary
    # ========================================================

    print()
    print("====================================================")
    print("FINAL PRIORITISATION SUMMARY")
    print("====================================================")

    print(
        "Total prioritised variants:",
        len(results)
    )

    print()

    for tier in range(1, 6):

        count = sum(
            1
            for row in results
            if row["TIER"] == tier
        )

        print(
            f"Tier {tier}: {count}"
        )

    print()
    print(
        "Output file:",
        args.output
    )

    print("====================================================")
    print()


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":

    main()
