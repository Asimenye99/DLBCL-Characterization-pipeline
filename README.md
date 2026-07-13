# DLBCL-Characterization-pipeline
Repository for functional characterization of DLBCL in Malawi

## Overview

This repository contains a Snakemake-based bioinformatics pipeline developed for the in silico functional characterisation of pathogenic variants in cancer-associated genes in diffuse large B-cell lymphoma (DLBCL) among Malawian patients.

The workflow integrates whole-exome sequencing (WES) data processing, including quality control, read trimming, alignment to the GRCh38 reference genome, variant calling, annotation, and prioritisation of clinically and biologically relevant variants. Downstream analyses include extraction of DLBCL-associated driver gene variants, protein sequence analysis, and preparation of variants for structural and functional assessment.

The pipeline is designed to provide a reproducible framework for identifying potentially pathogenic variants and investigating their possible effects on protein structure, stability, and disease-related molecular pathways.

```text
DLBCL-Characterization-pipeline/
│
├── workflow/
│   └── Snakefile                    # Main Snakemake workflow
│
├── config/
│   └── config.yaml                  # Pipeline configuration and parameters
│
├── scripts/
│   ├── extract_protein_sequences.py # Extract protein sequences for modelling
│   ├── mutate_sequence.py           # Introduce amino acid substitutions
│   └── prioritize_variants.py       # Variant filtering and prioritisation
│
├── reference/
│   └── adapters/
│       └── NexteraPE-PE.fa          # Adapter sequences for read trimming
│
├── data/
│   └── (local sequencing data - not included)
│
├── results/
│   └── (generated analysis outputs - not included)
│
├── logs/
│   └── (workflow execution logs - not included)
│
├── benchmarks/
│   └── (runtime and resource reports)
│
├── annotation/
│   └── (annotation resources - not included)
│
├── samples.tsv                      # Sample metadata and FASTQ paths
├── driver_genes.txt                 # DLBCL driver gene list
├── README.md                         # Documentation
└── .gitignore                        # Files excluded from version control
