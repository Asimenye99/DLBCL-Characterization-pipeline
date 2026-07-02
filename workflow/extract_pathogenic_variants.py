#!/usr/bin/env python3

import argparse
import csv
import re
from cyvcf2 import VCF
import pandas as pd

PATHOGENIC_RE = re.compile(r"pathogenic", re.IGNORECASE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="SAMPLE")
    args = ap.parse_args()

    vcf = VCF(args.vcf)

    # ----------------------------
    # Parse CSQ header from VEP
    # ----------------------------
    csq_header = None
    for h in vcf.header_iter():
        info = h.info(extra=True)
        if info.get("ID") == "CSQ":
            desc = info["Description"]
            csq_header = desc.split("Format:")[1].strip().strip('"').split("|")
            break

    if not csq_header:
        raise ValueError("CSQ header not found in VCF")

    idx = {k: i for i, k in enumerate(csq_header)}

    raw_rows = []

    # ----------------------------
    # Extract pathogenic annotations
    # ----------------------------
    for var in vcf:
        csq = var.INFO.get("CSQ")
        if not csq:
            continue

        for ann in csq.split(","):
            f = ann.split("|")

            if len(f) != len(csq_header):
                continue

            gene = f[idx["SYMBOL"]] if "SYMBOL" in idx else ""
            clin = f[idx["CLIN_SIG"]] if "CLIN_SIG" in idx else ""
            impact = f[idx["IMPACT"]] if "IMPACT" in idx else ""
            hgvsp = f[idx["HGVSp"]] if "HGVSp" in idx else ""

            if not clin:
                continue

            if not PATHOGENIC_RE.search(clin):
                continue

            raw_rows.append([
                args.sample,
                var.CHROM,
                var.POS,
                var.REF,
                ",".join(var.ALT),
                gene,
                impact,
                clin,
                hgvsp
            ])

    if not raw_rows:
        print("No pathogenic variants found")
        return

    # ----------------------------
    # Save RAW transcript-level results
    # ----------------------------
    raw_file = args.out.replace(".tsv", ".raw.tsv")

    with open(raw_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["SAMPLE","CHROM","POS","REF","ALT","GENE","IMPACT","CLIN_SIG","HGVSp"])
        w.writerows(raw_rows)

    # ----------------------------
    # Collapse duplicates → UNIQUE VARIANTS
    # ----------------------------
    df = pd.DataFrame(raw_rows, columns=[
        "SAMPLE","CHROM","POS","REF","ALT","GENE","IMPACT","CLIN_SIG","HGVSp"
    ])

    unique = df.drop_duplicates(subset=["CHROM","POS","REF","ALT"])

    unique_file = args.out.replace(".tsv", ".unique.tsv")
    unique.to_csv(unique_file, sep="\t", index=False)

    print(f"Transcript-level pathogenic annotations: {len(df)}")
    print(f"Unique pathogenic variants: {len(unique)}")
    print(f"Saved raw → {raw_file}")
    print(f"Saved unique → {unique_file}")


if __name__ == "__main__":
    main()
