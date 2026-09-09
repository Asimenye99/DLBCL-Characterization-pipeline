#!/usr/bin/env python3

import argparse
import requests
import pandas as pd
import sys

GP_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def main():
    parser = argparse.ArgumentParser(
        description="KEGG pathway enrichment using g:Profiler."
    )

    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    try:
        # Read gene list
        with open(args.genes) as f:
            genes = [x.strip() for x in f if x.strip()]

        if not genes:
            raise ValueError("Gene list is empty.")

        print(f"Genes submitted: {len(genes)}")

        # g:Profiler query
        payload = {
            "organism": "hsapiens",
            "query": genes,
            "sources": ["KEGG"],
            "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS"
        }

        response = requests.post(
            GP_URL,
            json=payload,
            timeout=300
        )

        response.raise_for_status()

        data = response.json()
        results = data.get("result", [])

        if not results:
            print("No significant KEGG pathways found.")

            pd.DataFrame().to_csv(
                args.output,
                sep="\t",
                index=False
            )

            return

        rows = []

        for r in results:

            # g:Profiler may return intersections as
            # strings or nested lists.
            intersections = r.get("intersections", [])

            formatted_intersections = ",".join(
                item[0]
                if isinstance(item, list) and item
                else str(item)
                for item in intersections
            )

            rows.append({
                "term_id": r.get("native"),
                "pathway": r.get("name"),
                "p_value": r.get("p_value"),
                "significant": r.get("significant"),
                "intersection_size": r.get("intersection_size"),
                "term_size": r.get("term_size"),
                "query_size": r.get("query_size"),
                "intersection": formatted_intersections
            })

        df = pd.DataFrame(rows).sort_values("p_value")

        df.to_csv(
            args.output,
            sep="\t",
            index=False
        )

        print(f"KEGG pathways returned: {len(df)}")
        print(f"Output: {args.output}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
