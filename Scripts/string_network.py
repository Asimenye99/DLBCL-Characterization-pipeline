#!/usr/bin/env python3
"""Retrieve the STRING protein-protein interaction network for a gene list.

The identifiers are sent in the body of a POST request, separated by carriage
returns, as documented for the STRING API. (Passing them in the URL query
string, pre-encoded as %0d, is double-encoded by `requests` to %250d, which
STRING cannot parse and which also fails for long gene lists.)
"""
import argparse
import sys
import time

import requests

STRING_URL = "https://string-db.org/api/tsv/network"


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve STRING protein-protein interaction network."
    )
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--context-genes", default=None,
        help="Optional gene list (e.g. the Reddy DLBCL drivers) submitted "
             "together with the candidates so that candidate-to-driver "
             "associations can be examined.",
    )
    parser.add_argument(
        "--context-output", default=None,
        help="Edges with at least one candidate endpoint, from the combined query.",
    )
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--required-score", type=int, default=400)
    args = parser.parse_args()

    try:
        with open(args.genes) as f:
            genes = [x.strip() for x in f if x.strip()]
        if not genes:
            raise ValueError("Gene list is empty.")

        # Remove duplicates, preserving order
        genes = list(dict.fromkeys(genes))
        print(f"Genes submitted: {len(genes)}")

        candidates = set(genes)
        query = list(genes)
        if args.context_genes:
            with open(args.context_genes) as f:
                query += [x.strip() for x in f if x.strip()]
            query = list(dict.fromkeys(query))
            print(f"Genes submitted incl. context: {len(query)}")

        data = {
            "identifiers": "\r".join(query),
            "species": args.species,
            "caller_identity": "DLBCL_Analysis",
            "required_score": args.required_score,
        }

        print(f"Querying STRING with {len(genes)} genes...")

        response = None
        for attempt in range(1, 4):
            try:
                response = requests.post(STRING_URL, data=data, timeout=300)
                response.raise_for_status()
                break
            except requests.RequestException as err:
                print(f"Attempt {attempt} failed: {err}", file=sys.stderr)
                if attempt == 3:
                    raise
                time.sleep(10 * attempt)

        full = response.text
        if not full.strip():
            raise ValueError("STRING returned an empty response.")

        # Candidate-only network (both endpoints are candidate genes) keeps
        # --output comparable with the original analysis.
        flines = full.strip().splitlines()
        header = flines[0].split("\t")
        ia, ib = header.index("preferredName_A"), header.index("preferredName_B")
        rows = [l for l in flines[1:] if l.split("\t")[ia] in candidates
                or l.split("\t")[ib] in candidates]
        both = [l for l in rows if l.split("\t")[ia] in candidates
                and l.split("\t")[ib] in candidates]
        if args.context_output:
            with open(args.context_output, "w") as f:
                f.write("\n".join([flines[0]] + rows) + "\n")
        text = "\n".join([flines[0]] + both) + "\n"

        lines = text.strip().splitlines()
        print(f"STRING returned {len(lines) - 1} interactions "
              f"(score >= {args.required_score / 1000:.3f}).")

        with open(args.output, "w") as f:
            f.write(text)
        print(f"STRING network written to: {args.output}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
