#!/usr/bin/env python3

import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Draw STRING protein-protein interaction network."
    )
    parser.add_argument("--network", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.network, sep="\t")

        if df.empty:
            raise ValueError("STRING network is empty.")

        print(f"Interactions: {len(df)}")

        G = nx.Graph()

        for _, row in df.iterrows():

            protein_a = row.get("preferredName_A")
            protein_b = row.get("preferredName_B")
            score = row.get("score", 0)

            if pd.notna(protein_a) and pd.notna(protein_b):
                G.add_edge(
                    protein_a,
                    protein_b,
                    weight=float(score)
                )

        if len(G.nodes) == 0:
            raise ValueError("No valid network nodes found.")

        plt.figure(figsize=(12, 10))

        pos = nx.spring_layout(
            G,
            seed=42,
            k=1.5
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            node_size=1000
        )

        nx.draw_networkx_edges(
            G,
            pos,
            alpha=0.5
        )

        nx.draw_networkx_labels(
            G,
            pos,
            font_size=9
        )

        plt.title("DLBCL Candidate Gene Protein-Protein Interaction Network")

        plt.axis("off")
        plt.tight_layout()

        plt.savefig(
            args.output,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        print(f"Network figure: {args.output}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
