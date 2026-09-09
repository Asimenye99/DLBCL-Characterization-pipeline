#!/usr/bin/env python3

import argparse
import pandas as pd
import base64
import sys
from pathlib import Path


def table_html(path, title):
    try:
        df = pd.read_csv(path, sep="\t")

        if df.empty:
            return f"<h2>{title}</h2><p>No significant results.</p>"

        return (
            f"<h2>{title}</h2>"
            + df.head(20).to_html(
                index=False,
                classes="results"
            )
        )

    except Exception as e:
        return f"<h2>{title}</h2><p>Error: {e}</p>"


def image_base64(path):
    with open(path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode()

    return encoded


def main():
    parser = argparse.ArgumentParser(
        description="Build DLBCL pathway analysis HTML report."
    )

    parser.add_argument("--network", required=True)
    parser.add_argument("--reactome", required=True)
    parser.add_argument("--kegg", required=True)
    parser.add_argument("--go", required=True)
    parser.add_argument("--figure", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    try:
        network_df = pd.read_csv(args.network, sep="\t")

        image = image_base64(args.figure)

        html = f"""
<!DOCTYPE html>

<html>

<head>

<title>DLBCL Pathway Analysis</title>

<style>

body {{
    font-family: Arial, sans-serif;
    margin: 40px;
}}

h1 {{
    margin-bottom: 30px;
}}

h2 {{
    margin-top: 40px;
}}

.results {{
    border-collapse: collapse;
    width: 100%;
}}

.results th,
.results td {{
    border: 1px solid #ddd;
    padding: 6px;
    font-size: 12px;
}}

.results th {{
    background-color: #f2f2f2;
}}

img {{
    max-width: 100%;
}}

</style>

</head>

<body>

<h1>DLBCL Candidate Gene Pathway Analysis</h1>

<h2>Protein-Protein Interaction Network</h2>

<p>
Number of STRING interactions:
{len(network_df)}
</p>

<img src="data:image/png;base64,{image}">

{table_html(args.go, "GO Enrichment")}

{table_html(args.kegg, "KEGG Pathway Enrichment")}

{table_html(args.reactome, "Reactome Pathway Enrichment")}

</body>

</html>
"""

        Path(args.output).write_text(html)

        print(f"Report written to: {args.output}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
