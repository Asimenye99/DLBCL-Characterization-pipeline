#!/usr/bin/env python3
"""Network context for variant-bearing genes (STRING).

STRING holds known and predicted associations. It does not show that a variant
has changed an interaction. This script therefore describes only:
  1. connectivity of each candidate gene in the candidate network (degree,
     weighted degree, betweenness; hub = top decile of degree),
  2. associations of each candidate with Reddy et al. (2017) DLBCL drivers,
  3. the FoldX classification of the candidate's missense variants, so that
     structurally affected proteins can be located in the network.
"""
import argparse
import networkx as nx
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--network", required=True, help="candidate-only STRING edges")
ap.add_argument("--context-network", required=True, help="edges incl. Reddy genes")
ap.add_argument("--candidate-genes", required=True, help="candidate_genes.tsv")
ap.add_argument("--reddy-genes", required=True, help="comma-separated")
ap.add_argument("--structural", default=None, help="DLBCL_structural_effects.tsv")
ap.add_argument("--output", required=True)
a = ap.parse_args()

cand = pd.read_csv(a.candidate_genes, sep="\t", dtype=str)
genes = set(cand.gene)
reddy = {g.strip() for g in a.reddy_genes.split(",") if g.strip()}

def graph(path):
    e = pd.read_csv(path, sep="\t")
    g = nx.Graph()
    for r in e.itertuples():
        g.add_edge(r.preferredName_A, r.preferredName_B, weight=float(r.score))
    return g

G = graph(a.network)
C = graph(a.context_network)
deg = dict(G.degree())
wdeg = dict(G.degree(weight="weight"))
btw = nx.betweenness_centrality(G) if len(G) else {}
cut = pd.Series(list(deg.values())).quantile(0.9) if deg else 0

rows = []
for g in sorted(genes):
    nb = set(C.neighbors(g)) if g in C else set()
    rn = sorted(nb & reddy)
    rows.append({"gene": g,
                 "GENE_SET": cand.loc[cand.gene == g, "GENE_SET"].iloc[0],
                 "degree_candidate_network": deg.get(g, 0),
                 "weighted_degree": round(wdeg.get(g, 0.0), 3),
                 "betweenness": round(btw.get(g, 0.0), 4),
                 "hub_top_decile": deg.get(g, 0) >= cut and deg.get(g, 0) > 0,
                 "n_reddy_neighbours": len(rn),
                 "reddy_neighbours": ",".join(rn)})
out = pd.DataFrame(rows)

if a.structural:
    s = pd.read_csv(a.structural, sep="\t")
    gc = next((c for c in s.columns if c.lower() in ("gene", "gene_symbol")), None)
    cc = next((c for c in s.columns if "class" in c.lower() or "effect" in c.lower()), None)
    if gc and cc:
        m = s.groupby(gc)[cc].apply(lambda x: ",".join(sorted(set(map(str, x))))).rename("foldx_classes")
        out = out.merge(m, left_on="gene", right_index=True, how="left")
out.sort_values(["degree_candidate_network", "betweenness"], ascending=False).to_csv(a.output, sep="\t", index=False)
print(f"Candidates: {len(genes)}; in network: {sum(1 for g in genes if g in G)}; "
      f"hubs (degree >= {cut:.0f}): {int(out.hub_top_decile.sum())}; "
      f"with >=1 Reddy neighbour: {int((out.n_reddy_neighbours > 0).sum())}")
