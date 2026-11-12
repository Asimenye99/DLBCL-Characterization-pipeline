#!/usr/bin/env python3
"""
Build the gene lists and ID maps needed for structural modelling directly from
the prioritised variant table, so that no hand-curated gene panel is required.

Input : TSV from prioritize_variants.py (needs GENE, CONSEQUENCE, ENSP,
        SWISSPROT, CANONICAL; use the pathogenic-evidence table).
Output: --genes        text file, one gene per line: genes with at least one
                       missense variant and a usable UniProt + ENSP mapping
        --uniprot-map  TSV gene<TAB>uniprot_id   (SWISSPROT from VEP, version stripped)
        --ensp-map     TSV gene<TAB>ensp         (canonical transcript's protein)
        --unmapped     TSV of genes with missense variants that could not be
                       mapped (reported, not silently dropped)

The UniProt accession comes from VEP's SWISSPROT field (requires --uniprot,
included in --everything). Check that the ENSP and the AlphaFold canonical
isoform describe the same protein before trusting a structure; the
mutation-list step should confirm the wild-type residue at each position.
"""
import argparse
import csv
from collections import Counter, defaultdict


def first_accession(value):
    """'P04637.263&Q9XXXX.1' -> 'P04637' (reviewed Swiss-Prot entry, no version)."""
    for item in (value or "").replace(",", "&").split("&"):
        item = item.strip().split(".")[0]
        if item:
            return item
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--uniprot-map", required=True)
    ap.add_argument("--ensp-map", required=True)
    ap.add_argument("--unmapped", required=True)
    a = ap.parse_args()

    ensp_votes = defaultdict(Counter)    # gene -> Counter of ENSP (canonical weighted)
    uni_votes = defaultdict(Counter)
    missense_genes = set()

    with open(a.variants) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["GENE"].strip()
            if not gene or "missense_variant" not in row["CONSEQUENCE"].split("&"):
                continue
            missense_genes.add(gene)
            weight = 3 if row.get("CANONICAL") == "YES" else 1
            if row.get("ENSP"):
                ensp_votes[gene][row["ENSP"]] += weight
            acc = first_accession(row.get("SWISSPROT"))
            if acc:
                uni_votes[gene][acc] += weight

    mapped, unmapped = [], []
    for gene in sorted(missense_genes):
        if ensp_votes[gene] and uni_votes[gene]:
            mapped.append((gene, uni_votes[gene].most_common(1)[0][0],
                           ensp_votes[gene].most_common(1)[0][0]))
        else:
            reason = ("no_ENSP" if not ensp_votes[gene] else "") + \
                     ("no_SWISSPROT" if not uni_votes[gene] else "")
            unmapped.append((gene, reason))

    with open(a.genes, "w") as f:
        f.write("".join(g + "\n" for g, _, _ in mapped))
    with open(a.uniprot_map, "w") as f:
        f.write("gene\tuniprot_id\n" + "".join(f"{g}\t{u}\n" for g, u, _ in mapped))
    with open(a.ensp_map, "w") as f:
        f.write("gene\tensp\n" + "".join(f"{g}\t{e}\n" for g, _, e in mapped))
    with open(a.unmapped, "w") as f:
        f.write("gene\treason\n" + "".join(f"{g}\t{r}\n" for g, r in unmapped))
    print(f"Genes with missense variants: {len(missense_genes)}; mapped: {len(mapped)}; unmapped: {len(unmapped)}")


if __name__ == "__main__":
    main()
