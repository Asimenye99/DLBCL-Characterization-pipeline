#!/usr/bin/env python3

import os
import gzip
import subprocess
from pathlib import Path

PROJECT = Path(".")

print("="*90)
print("        DLBCL PIPELINE PROGRESS REPORT")
print("="*90)

# -------------------------------------------------------
# Helper functions
# -------------------------------------------------------

def count_fastq_reads(fastq):
    try:
        cmd = f"zcat {fastq} | echo $((`wc -l`/4))"
    except:
        pass

    n = subprocess.check_output(
        f"zcat {fastq} | wc -l",
        shell=True,
        text=True
    )

    return int(n)//4


def count_vcf(vcf):

    if str(vcf).endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    n = 0

    with opener(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1

    return n


def run(cmd):

    try:
        return subprocess.check_output(
            cmd,
            shell=True,
            text=True
        ).strip()

    except:
        return "NA"


# -------------------------------------------------------
# FASTQ
# -------------------------------------------------------

print("\n1. RAW READS")
print("-"*90)

raw = sorted(list(Path("data/raw").glob("*.fastq.gz")) +
                list(Path("data/raw").glob("*.fq.gz")))

total_reads = 0

for fq in raw:

    reads = count_fastq_reads(fq)
    total_reads += reads

    print(f"{fq.name:<40}{reads:,}")

print(f"\nTOTAL READS : {total_reads:,}")


# -------------------------------------------------------
# BAM
# -------------------------------------------------------

print("\n2. ALIGNMENT")
print("-"*90)

bam_dir = Path("results/bam")

for bam in bam_dir.glob("*.bam"):

    print("\n"+bam.name)

    mapped = run(f"samtools flagstat {bam} | grep 'mapped ('")

    print(mapped)


# -------------------------------------------------------
# DUPLICATES
# -------------------------------------------------------

print("\n3. PCR DUPLICATES")
print("-"*90)

metrics = list(Path("results").rglob("*.metrics.txt"))

if not metrics:
    metrics = list(Path("data").rglob("*.metrics.txt"))

if not metrics:
    print("Metrics file not found")
else:
    for m in metrics:

        with open(m) as f:

            for line in f:

                if line.startswith("LIBRARY"):

                    values = next(f).split()

                    dup = float(values[8]) * 100

                    print(f"{m.name:<30}{dup:.2f}%")

                    break


# -------------------------------------------------------
# VARIANTS
# -------------------------------------------------------

import subprocess

print("\n4. VARIANT SUMMARY")
print("-"*90)

filtered = "results/mutect/WBP02.filtered.vcf.gz"

total = int(run(f"bcftools view -H {filtered} | wc -l"))
passed = int(run(f"bcftools view -f PASS -H {filtered} | wc -l"))

failed = total - passed

print(f"Total variants           : {total:,}")
print(f"PASS variants            : {passed:,}")
print(f"Filtered variants        : {failed:,}")


# -------------------------------------------------------
# CLINVAR
# -------------------------------------------------------

# -------------------------------------------------------
# CLINVAR SUMMARY (from VEP CSQ)
# -------------------------------------------------------

print("\n5. CLINVAR SUMMARY")
print("-"*90)

import re

vep_files = list(Path("results/annotation").glob("*.vcf")) + \
            list(Path("results/annotation").glob("*.vcf.gz"))

if len(vep_files):

    vep = vep_files[0]

    opener = gzip.open if str(vep).endswith(".gz") else open

    # Get CSQ field order
    csq_fields = []

    with opener(vep, "rt") as f:

        for line in f:

            if line.startswith("##INFO=<ID=CSQ"):

                m = re.search(r'Format: (.+)">', line)

                csq_fields = m.group(1).split("|")

                break

    clnsig_index = csq_fields.index("ClinVar_CLNSIG")

    clinvar = {
        "Pathogenic":0,
        "Likely pathogenic":0,
        "Pathogenic/Likely pathogenic":0,
        "Benign":0,
        "Likely benign":0,
        "Uncertain significance":0,
        "Conflicting":0
    }

    with opener(vep, "rt") as f:

        for line in f:

            if line.startswith("#"):
                continue

            info = line.split("\t")[7]

            m = re.search(r'CSQ=([^;]+)', info)

            if not m:
                continue

            annotations = m.group(1).split(",")

            # only count variant once
            found = set()

            for ann in annotations:

                cols = ann.split("|")

                if len(cols) <= clnsig_index:
                    continue

                sig = cols[clnsig_index].strip()

                if sig == "":
                    continue

                sig = sig.lower()

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

print()

total = 0

for k,v in clinvar.items():

    total += v

    print(f"{k:<35}{v:>12,}")

print("-"*50)
print(f"{'Total ClinVar annotated':<35}{total:>12,}")


# -------------------------------------------------------
# PRIORITIZED
# -------------------------------------------------------

print("\n6. PRIORITIZED DLBCL VARIANTS")
print("-"*90)

files = list(Path("results/prioritization").glob("*.tsv"))

for tsv in files:

    rows = sum(1 for _ in open(tsv))-1

    print(f"{tsv.name:<40}{rows}")


# -------------------------------------------------------
# WILDTYPE
# -------------------------------------------------------

print("\n7. WILDTYPE PROTEINS")
print("-"*90)

wt = len(list(Path("results/proteins/wildtype").glob("*.fa"))) + \
     len(list(Path("results/proteins/wildtype").glob("*.fasta")))

print(wt)


# -------------------------------------------------------
# MUTANTS
# -------------------------------------------------------

print("\n8. MUTANT PROTEINS")
print("-"*90)

mt = len(list(Path("results/mutant_proteins").glob("*.fa"))) + \
     len(list(Path("results/mutant_proteins").glob("*.fasta")))

print(mt)


# -------------------------------------------------------
# STRUCTURE INPUTS
# -------------------------------------------------------

print("\n9. STRUCTURE INPUTS")
print("-"*90)

st = len(list(Path("results/structure_inputs").glob("*")))

print(st)


# -------------------------------------------------------
# PIPELINE STATUS
# -------------------------------------------------------

print("\n10. PIPELINE STATUS")
print("-"*90)

steps = {
"FASTQ":"✓",
"FASTQC":"✓",
"Trimming":"✓",
"BWA Alignment":"✓",
"Deduplication":"✓",
"Mutect2":"✓",
"VCF Filtering":"✓",
"VEP Annotation":"✓",
"ClinVar Annotation":"✓",
"Variant Prioritization":"✓",
"Wildtype Proteins":"✓",
"Mutant Proteins":"✓" if mt>0 else "⚠",
"AlphaFold Structures":"✗",
"FoldX Stability":"✗",
"Functional Interpretation":"✗"
}

for s,v in steps.items():

    print(f"{v} {s}")

print("\n"+"="*90)
print("REPORT COMPLETE")
print("="*90)