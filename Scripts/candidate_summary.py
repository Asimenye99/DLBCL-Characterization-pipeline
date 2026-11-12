#!/usr/bin/env python3
"""Funnel counts, gene-set overlap and the candidate variant table.

Inputs  : the all-variants table (DLBCL_candidate_variants.tsv) and the
          pathogenic-evidence table, both written by prioritize_variants.py
          with --census, plus the Census and Reddy gene lists.
Outputs : funnel_summary.tsv, gene_set_overlap.tsv, candidate_variant_table.tsv
Gene membership is an annotation; it is never used as evidence of pathogenicity.
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--all-variants", required=True)
ap.add_argument("--pathogenic", required=True)
ap.add_argument("--census", required=True)
ap.add_argument("--reddy-genes", required=True, help="comma-separated symbols")
ap.add_argument("--outdir", required=True)
ap.add_argument("--levels", nargs="+", default=["1", "2", "3"], help="evidence levels in the candidate set")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

allv = pd.read_csv(a.all_variants, sep="\t", dtype=str, keep_default_na=False)
path = pd.read_csv(a.pathogenic, sep="\t", dtype=str, keep_default_na=False)
cen = pd.read_csv(a.census, dtype=str, keep_default_na=False)
cgc = set(cen["Gene Symbol"])
reddy = {g.strip() for g in a.reddy_genes.split(",") if g.strip()}

for d in (allv, path):
    d["VAR"] = d.CHROM + ":" + d.POS + d.REF + ">" + d.ALT
cand = path[(path.GENE_SET != "Neither") & path.EVIDENCE_LEVEL.isin(a.levels)].copy()
trunc_set = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
             "splice_donor_variant", "start_lost", "stop_lost"}
cand["CLASS"] = cand.CONSEQUENCE.apply(
    lambda c: "truncating" if trunc_set & set(c.split("&"))
    else ("missense" if "missense_variant" in c.split("&") else "other_protein_altering"))

def line(step, d):
    return {"step": step, "calls": len(d), "unique_variants": d.VAR.nunique(),
            "genes": d.GENE.nunique(), "patients": d.SAMPLE.nunique()}

rows = [line("Protein-altering calls retained after frequency filter", allv),
        line("With pathogenic evidence (Levels 1-3), any gene", path),
        line("ClinVar pathogenic/likely pathogenic (Level 1), any gene", path[path.EVIDENCE_LEVEL == "1"]),
        line(f"... in a Census or Reddy gene, level(s) {'+'.join(a.levels)} (candidate set)", cand)]
for gs in ("Both", "Census_only", "Reddy_only"):
    rows.append(line(f"      GENE_SET = {gs}", cand[cand.GENE_SET == gs]))
for lv in ("1", "2", "3"):
    rows.append(line(f"      Evidence Level {lv}", cand[cand.EVIDENCE_LEVEL == lv]))
for c in ("missense", "truncating", "other_protein_altering"):
    rows.append(line(f"      Class {c}", cand[cand.CLASS == c]))
for t in ("1", "2"):
    rows.append(line(f"      Census Tier {t}", cand[cand.CGC_TIER == t]))
pd.DataFrame(rows).to_csv(f"{a.outdir}/funnel_summary.tsv", sep="\t", index=False)

# gene-level overlap of the two reference lists, and which genes carry candidates
cg = set(cand.GENE)
ov = [("Census genes (list)", len(cgc)), ("Reddy genes (list, symbols as configured)", len(reddy)),
      ("Reddy genes also in Census (exact symbol)", len(reddy & cgc)),
      ("Candidate genes", len(cg)),
      ("  Census and Reddy", len(cg & cgc & reddy)),
      ("  Census only", len((cg & cgc) - reddy)),
      ("  Reddy only", len((cg & reddy) - cgc))]
pd.DataFrame(ov, columns=["item", "n"]).to_csv(f"{a.outdir}/gene_set_overlap.tsv", sep="\t", index=False)
pd.DataFrame({"gene": sorted(cg),
              "GENE_SET": [cand[cand.GENE == g].GENE_SET.iloc[0] for g in sorted(cg)],
              "CGC_TIER": [cand[cand.GENE == g].CGC_TIER.iloc[0] for g in sorted(cg)],
              "n_variants": [cand[cand.GENE == g].VAR.nunique() for g in sorted(cg)],
              "n_patients": [cand[cand.GENE == g].SAMPLE.nunique() for g in sorted(cg)]}
             ).to_csv(f"{a.outdir}/candidate_genes.tsv", sep="\t", index=False)

cols = ["GENE", "CHROM", "POS", "REF", "ALT", "ENSP", "HGVSc", "HGVSp", "CONSEQUENCE", "CLASS",
        "MAX_AF", "EVIDENCE_LEVEL", "PATHOGENIC_EVIDENCE", "ClinVar_CLNSIG", "CADD_PHRED", "REVEL",
        "AlphaMissense_score", "IN_CGC", "CGC_TIER", "IN_DLBCL_LIST", "GENE_SET", "SOURCE", "SAMPLE"]
agg = (cand.groupby(["VAR", "GENE", "HGVSp"], sort=False)
       .agg({**{c: "first" for c in cols if c not in ("SAMPLE", "SOURCE")},
             "SOURCE": lambda x: ";".join(sorted(set(";".join(x).split(";")))),
             "SAMPLE": lambda x: ",".join(sorted(set(x)))}).reset_index(drop=True))
agg = agg.rename(columns={"SAMPLE": "PATIENTS"})
agg["N_PATIENTS"] = agg.PATIENTS.str.count(",") + 1
agg.sort_values(["EVIDENCE_LEVEL", "GENE"]).to_csv(f"{a.outdir}/candidate_variant_table.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(ov, columns=["item", "n"]).to_string(index=False))
