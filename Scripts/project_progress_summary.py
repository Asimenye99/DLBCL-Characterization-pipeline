#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import gzip
import re
import subprocess
import sys


# ============================================================
# PROJECT CONFIGURATION
# ============================================================

PROJECT = Path("/head/NFS/dndovi/DLBCL_Analysis")

RESULTS = PROJECT / "results"

MUTECT_DIR = RESULTS / "mutect"
VEP_DIR = RESULTS / "annotation"
PRIORITY_DIR = RESULTS / "prioritization"
STRUCTURAL_FILE = RESULTS / "structural_effects" / "DLBCL_structural_effects.tsv"
PATHWAY_DIR = RESULTS / "pathways"

PROTEIN_DIR = RESULTS / "proteins"
MUTANT_DIR = RESULTS / "mutant_proteins"
STRUCTURE_DIR = RESULTS / "structures"
FOLDX_DIR = RESULTS / "foldx"

BAM_DIR = RESULTS / "bam"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def count_vcf(vcf):
    """Count variant records in a VCF/VCF.GZ without loading it into memory."""

    opener = gzip.open if str(vcf).endswith(".gz") else open

    count = 0

    try:
        with opener(vcf, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    count += 1
    except Exception:
        return 0

    return count


def count_tsv_rows(tsv):
    """Count data rows in TSV excluding header."""

    try:
        with open(tsv) as f:
            return max(0, sum(1 for _ in f) - 1)
    except Exception:
        return 0


def run(cmd):
    """Run shell command safely."""

    try:
        return subprocess.check_output(
            cmd,
            shell=True,
            text=True,
            stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return ""


def find_first(paths):
    """Return first existing file from list."""

    for path in paths:
        if path.exists():
            return path

    return None


# ============================================================
# HEADER
# ============================================================

print("=" * 100)
print("              DLBCL PIPELINE → THESIS RESULTS REPORT")
print("=" * 100)

print()
print("Project directory:")
print(f"  {PROJECT}")

print()


# ============================================================
# 1. COHORT / SAMPLE SUMMARY
# ============================================================

print("1. COHORT SUMMARY")
print("-" * 100)

bam_files = sorted(BAM_DIR.glob("*.marked.bam"))

normalized_vcfs = sorted(
    MUTECT_DIR.glob("*.normalized.vcf.gz")
)

vep_files = sorted(
    VEP_DIR.glob("*.vep.vcf")
)

candidate_files = sorted(
    PRIORITY_DIR.glob("*.tsv")
)

print(f"Samples with marked BAM files        : {len(bam_files)}")
print(f"Samples with normalized VCFs         : {len(normalized_vcfs)}")
print(f"Samples with VEP annotation          : {len(vep_files)}")

print()


# ============================================================
# 2. VARIANT CALLING
# ============================================================

print("2. VARIANT CALLING")
print("-" * 100)

filtered_vcfs = sorted(
    list(MUTECT_DIR.glob("*.filtered.vcf.gz")) +
    list(MUTECT_DIR.glob("*.filtered.vcf"))
)

total_variants = 0
total_pass = 0

if filtered_vcfs:

    for vcf in filtered_vcfs:

        total = count_vcf(vcf)

        passed = run(
            f"bcftools view -f PASS -H '{vcf}' | wc -l"
        )

        try:
            passed = int(passed)
        except Exception:
            passed = 0

        total_variants += total
        total_pass += passed

    print(f"Filtered Mutect2 VCFs               : {len(filtered_vcfs)}")
    print(f"Total variant records               : {total_variants:,}")
    print(f"PASS variant records                : {total_pass:,}")

else:

    print("No filtered Mutect2 VCFs found.")

print()


# ============================================================
# 3. VEP ANNOTATION
# ============================================================

print("3. FUNCTIONAL ANNOTATION — VEP")
print("-" * 100)

if vep_files:

    total_vep_variants = sum(
        count_vcf(vcf)
        for vcf in vep_files
    )

    print(f"VEP-annotated samples               : {len(vep_files)}")
    print(f"Annotated variant records            : {total_vep_variants:,}")

else:

    print("No VEP files found.")

print()


# ============================================================
# 4. CLINVAR
# ============================================================

print("4. CLINVAR CLINICAL SIGNIFICANCE")
print("-" * 100)

clinvar = Counter()

if vep_files:

    # Use all VEP files rather than only the first one
    for vep in vep_files:

        opener = gzip.open if str(vep).endswith(".gz") else open

        csq_fields = []

        try:

            with opener(vep, "rt") as f:

                for line in f:

                    if line.startswith("##INFO=<ID=CSQ"):

                        m = re.search(
                            r"Format: (.+)\">",
                            line
                        )

                        if m:
                            csq_fields = m.group(1).split("|")

                        break

            if "ClinVar_CLNSIG" not in csq_fields:
                continue

            clnsig_index = csq_fields.index(
                "ClinVar_CLNSIG"
            )

            with opener(vep, "rt") as f:

                for line in f:

                    if line.startswith("#"):
                        continue

                    info = line.split("\t")[7]

                    m = re.search(
                        r"CSQ=([^;]+)",
                        info
                    )

                    if not m:
                        continue

                    found = set()

                    for ann in m.group(1).split(","):

                        cols = ann.split("|")

                        if len(cols) <= clnsig_index:
                            continue

                        sig = cols[clnsig_index].strip().lower()

                        if not sig:
                            continue

                        if "pathogenic/likely_pathogenic" in sig:
                            found.add("Pathogenic/Likely pathogenic")

                        elif sig == "pathogenic":
                            found.add("Pathogenic")

                        elif sig == "likely_pathogenic":
                            found.add("Likely pathogenic")

                        elif sig == "benign":
                            found.add("Benign")

                        elif sig == "likely_benign":
                            found.add("Likely benign")

                        elif "uncertain" in sig:
                            found.add("Uncertain significance")

                        elif "conflicting" in sig:
                            found.add("Conflicting")

                    for item in found:
                        clinvar[item] += 1

        except Exception:
            continue

    for key in [
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
        "Benign",
        "Likely benign",
        "Uncertain significance",
        "Conflicting"
    ]:

        if clinvar[key] > 0:
            print(
                f"{key:<35} {clinvar[key]:>10,}"
            )

    print(
        f"{'Total ClinVar classifications':<35}"
        f"{sum(clinvar.values()):>10,}"
    )

else:

    print("VEP annotation unavailable.")

print()


# ============================================================
# 5. PRIORITIZED CANDIDATE VARIANTS
# ============================================================

print("5. PRIORITIZED DLBCL CANDIDATE VARIANTS")
print("-" * 100)

candidate_file = None

if candidate_files:

    # Prefer the main DLBCL candidate file
    preferred = PRIORITY_DIR / "DLBCL_candidate_variants.tsv"

    if preferred.exists():
        candidate_file = preferred

    else:
        candidate_file = candidate_files[0]

    candidate_count = count_tsv_rows(candidate_file)

    print(f"Candidate variant file              : {candidate_file.name}")
    print(f"Prioritized candidate variants      : {candidate_count:,}")

else:

    print("No prioritized candidate file found.")

print()


# ============================================================
# 6. CANDIDATE GENE SUMMARY
# ============================================================

print("6. CANCER-ASSOCIATED GENE SUMMARY")
print("-" * 100)

gene_counts = Counter()

if candidate_file:

    try:

        with open(candidate_file) as f:

            header = f.readline().strip().split("\t")

            possible_gene_columns = [
                "SYMBOL",
                "Gene",
                "GENE",
                "gene",
                "Gene_Symbol"
            ]

            gene_col = None

            for name in possible_gene_columns:

                if name in header:
                    gene_col = header.index(name)
                    break

            if gene_col is not None:

                for line in f:

                    cols = line.rstrip("\n").split("\t")

                    if len(cols) > gene_col:

                        gene = cols[gene_col].strip()

                        if gene:
                            gene_counts[gene] += 1

    except Exception as e:

        print(f"Could not read candidate file: {e}")


if gene_counts:

    print(f"Unique candidate genes              : {len(gene_counts)}")

    print()
    print("Top candidate genes:")
    print()

    for gene, count in gene_counts.most_common(15):

        print(
            f"  {gene:<20} {count:>6} variants"
        )

else:

    print("No gene-level candidate information found.")

print()


# ============================================================
# 7. STRUCTURAL ANALYSIS
# ============================================================

print("7. STRUCTURAL STABILITY ANALYSIS — FoldX")
print("-" * 100)

structural_counts = Counter()
structural_records = []

if STRUCTURAL_FILE.exists():

    with open(STRUCTURAL_FILE) as f:

        header = f.readline().strip().split("\t")

        for line in f:

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 5:
                continue

            gene = cols[0]
            variant = cols[1]
            pdb = cols[2]

            try:
                ddg = float(cols[3])
            except Exception:
                continue

            classification = cols[4]

            structural_counts[classification] += 1

            structural_records.append(
                (
                    gene,
                    variant,
                    pdb,
                    ddg,
                    classification
                )
            )

    total_structural = len(structural_records)

    print(
        f"Structurally analysed variants      : "
        f"{total_structural:,}"
    )

    print()

    print(
        f"Destabilizing                        : "
        f"{structural_counts['destabilizing']:,}"
    )

    print(
        f"Neutral                              : "
        f"{structural_counts['neutral']:,}"
    )

    print(
        f"Stabilizing                          : "
        f"{structural_counts['stabilizing']:,}"
    )

else:

    print("Structural effects file not found.")

print()


# ============================================================
# 8. TOP DESTABILIZING MUTATIONS
# ============================================================

print("8. TOP STRUCTURALLY DESTABILIZING VARIANTS")
print("-" * 100)

destabilizing = [
    x for x in structural_records
    if x[4] == "destabilizing"
]

destabilizing.sort(
    key=lambda x: x[3],
    reverse=True
)

if destabilizing:

    print(
        f"{'Gene':<12}"
        f"{'Variant':<12}"
        f"{'ddG (kcal/mol)':>18}"
    )

    print("-" * 45)

    for gene, variant, pdb, ddg, classification in destabilizing[:15]:

        print(
            f"{gene:<12}"
            f"{variant:<12}"
            f"{ddg:>18.3f}"
        )

else:

    print("No destabilizing variants found.")

print()


# ============================================================
# 9. DESTABILIZING VARIANTS BY GENE
# ============================================================

print("9. GENES WITH STRUCTURALLY DESTABILIZING VARIANTS")
print("-" * 100)

destab_gene_counts = Counter(
    x[0] for x in destabilizing
)

if destab_gene_counts:

    print(
        f"{'Gene':<20}"
        f"{'Destabilizing variants':>25}"
    )

    print("-" * 50)

    for gene, count in destab_gene_counts.most_common():

        print(
            f"{gene:<20}"
            f"{count:>25}"
        )

else:

    print("No destabilizing genes identified.")

print()


# ============================================================
# 10. PROTEIN / STRUCTURE SUMMARY
# ============================================================

print("10. PROTEIN & STRUCTURE ANALYSIS")
print("-" * 100)

wt_files = (
    list(PROTEIN_DIR.rglob("*.fa")) +
    list(PROTEIN_DIR.rglob("*.fasta"))
)

mutant_files = (
    list(MUTANT_DIR.glob("*.fa")) +
    list(MUTANT_DIR.glob("*.fasta"))
)

pdb_files = list(STRUCTURE_DIR.glob("*.pdb"))

print(f"Wildtype protein sequences           : {len(wt_files)}")
print(f"Mutant protein sequences             : {len(mutant_files)}")
print(f"Protein structures                   : {len(pdb_files)}")

print()


# ============================================================
# 11. PATHWAY ANALYSIS
# ============================================================

print("11. PATHWAY ANALYSIS")
print("-" * 100)

reactome = PATHWAY_DIR / "reactome_enrichment.tsv"
kegg = PATHWAY_DIR / "kegg_enrichment.tsv"
go = PATHWAY_DIR / "go_enrichment.tsv"

pathway_results = {
    "GO enrichment": go,
    "KEGG enrichment": kegg,
    "Reactome enrichment": reactome
}

for name, path in pathway_results.items():

    if path.exists():

        rows = count_tsv_rows(path)

        print(
            f"{name:<25} {rows:>8} results"
        )

    else:

        print(
            f"{name:<25} NOT FOUND"
        )

print()


# ============================================================
# 12. PATHWAY TOP RESULTS
# ============================================================

print("12. TOP PATHWAY RESULTS")
print("-" * 100)


def show_top_pathways(path, name, n=10):

    if not path.exists():
        return

    print()
    print(name)
    print("-" * 60)

    try:

        with open(path) as f:

            header = f.readline().strip().split("\t")

            # Identify likely description/name and p-value columns
            print(
                "Columns: " +
                " | ".join(header[:8])
            )

            for i, line in enumerate(f):

                if i >= n:
                    break

                print("  " + line.rstrip())

    except Exception as e:

        print(f"Unable to read {name}: {e}")


show_top_pathways(
    reactome,
    "Reactome",
    5
)

show_top_pathways(
    kegg,
    "KEGG",
    5
)

show_top_pathways(
    go,
    "GO",
    5
)

print()


# ============================================================
# 13. KEY THESIS NUMBERS
# ============================================================

print("=" * 100)
print("                    THESIS-READY KEY NUMBERS")
print("=" * 100)

print()

print("COHORT")
print(
    f"  Samples processed                  : "
    f"{len(bam_files)}"
)

print()

print("VARIANTS")

if filtered_vcfs:

    print(
        f"  PASS variant records               : "
        f"{total_pass:,}"
    )

print()

if candidate_file:

    print(
        f"  Prioritized candidate variants     : "
        f"{count_tsv_rows(candidate_file):,}"
    )

print(
    f"  Candidate genes                    : "
    f"{len(gene_counts)}"
)

print()

print("STRUCTURAL EFFECTS")

if structural_records:

    print(
        f"  Variants structurally analysed     : "
        f"{len(structural_records)}"
    )

    print(
        f"  Destabilizing                      : "
        f"{structural_counts['destabilizing']}"
    )

    print(
        f"  Neutral                            : "
        f"{structural_counts['neutral']}"
    )

    print(
        f"  Stabilizing                        : "
        f"{structural_counts['stabilizing']}"
    )

print()

print("PATHWAYS")

for name, path in pathway_results.items():

    status = "AVAILABLE" if path.exists() else "PENDING"

    print(
        f"  {name:<25}: {status}"
    )

print()


# ============================================================
# 14. PIPELINE COMPLETION
# ============================================================

print("=" * 100)
print("                    PIPELINE STATUS")
print("=" * 100)

steps = {

    "Alignment":
        len(bam_files) > 0,

    "Mutect2 variant calling":
        len(normalized_vcfs) > 0,

    "VEP annotation":
        len(vep_files) > 0,

    "Variant prioritization":
        candidate_file is not None,

    "Protein sequence generation":
        len(wt_files) > 0,

    "Mutant protein generation":
        len(mutant_files) > 0,

    "Structural analysis":
        STRUCTURAL_FILE.exists()
        and len(structural_records) > 0,

    "GO enrichment":
        go.exists(),

    "KEGG enrichment":
        kegg.exists(),

    "Reactome enrichment":
        reactome.exists(),

    "Pathway summary":
        (PATHWAY_DIR / "pathway_summary.html").exists()
}


for step, status in steps.items():

    symbol = "✓" if status else "⚠"

    print(
        f"{symbol} {step}"
    )

print()


# ============================================================
# 15. THESIS ANALYSIS STATUS
# ============================================================

print("=" * 100)
print("                    THESIS ANALYSIS STATUS")
print("=" * 100)

thesis = {

    "Variant landscape":
        total_pass > 0,

    "Functional annotation":
        len(vep_files) > 0,

    "Clinical significance":
        sum(clinvar.values()) > 0,

    "Candidate variant prioritization":
        candidate_file is not None,

    "Candidate gene analysis":
        len(gene_counts) > 0,

    "Structural stability":
        len(structural_records) > 0,

    "GO enrichment":
        go.exists(),

    "KEGG enrichment":
        kegg.exists(),

    "Reactome enrichment":
        reactome.exists()
}


for section, status in thesis.items():

    print(
        f"{'READY' if status else 'PENDING':<10} "
        f"{section}"
    )

print()


# ============================================================
# FINAL
# ============================================================

print("=" * 100)
print("REPORT COMPLETE")
print("=" * 100)

print()
print("NEXT ANALYSIS FOCUS")
print()
print("1. Variant landscape")
print("2. Pathogenicity / clinical significance")
print("3. Candidate DLBCL genes")
print("4. Structurally destabilizing mutations")
print("5. Pathway enrichment")
print("6. Biological interpretation")
print()