#!/usr/bin/env python3
"""Background gene universe for enrichment: Census genes plus Reddy genes.
Candidates are drawn from these genes, so enrichment is tested against them
rather than against the whole genome (which returns generic cancer terms)."""
import argparse, csv
ap = argparse.ArgumentParser()
ap.add_argument("--census", required=True); ap.add_argument("--reddy-genes", required=True)
ap.add_argument("--output", required=True); ap.add_argument("--reddy-output", default=None); a = ap.parse_args()
genes = {r["Gene Symbol"].strip() for r in csv.DictReader(open(a.census, newline="")) if r["Gene Symbol"].strip()}
genes |= {g.strip() for g in a.reddy_genes.split(",") if g.strip()}
open(a.output, "w").write("\n".join(sorted(genes)) + "\n")
if a.reddy_output:
    open(a.reddy_output, "w").write("\n".join(sorted(g.strip() for g in a.reddy_genes.split(",") if g.strip())) + "\n")
print(len(genes), "background genes")
