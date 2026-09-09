#!/usr/bin/env python3

import argparse
import pandas as pd
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Extract unique gene symbols from DLBCL candidate variants."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input, sep="\t")

        print(f"Input rows: {len(df)}")
        print(f"Columns: {', '.join(df.columns)}")

        # Possible gene column names
        candidates = [
            "SYMBOL",
            "Gene",
            "GENE",
            "gene",
            "Gene_Symbol",
            "gene_symbol",
            "SYMBOLS",
            "HGNC"
        ]

        gene_col = None

        for col in candidates:
            if col in df.columns:
                gene_col = col
                break

        if gene_col is None:
            raise ValueError(
                "Could not find a gene column. "
                "Expected one of: " + ", ".join(candidates)
            )

        genes = (
            df[gene_col]
            .dropna()
            .astype(str)
            .str.strip()
        )

        # Remove empty values
        genes = genes[genes != ""]

        # Remove values that may contain multiple genes
        expanded = []

        for gene in genes:
            for g in gene.replace(",", ";").split(";"):
                g = g.strip()
                if g:
                    expanded.append(g)

        unique_genes = sorted(set(expanded))

        with open(args.output, "w") as f:
            for gene in unique_genes:
                f.write(gene + "\n")

        print(f"Gene column: {gene_col}")
        print(f"Unique genes: {len(unique_genes)}")
        print(f"Output: {args.output}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
