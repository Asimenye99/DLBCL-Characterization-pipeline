#!/usr/bin/env python3
"""
extract_protein_sequences.py

Fetches canonical wild-type protein sequences (FASTA) for each ENSP ID
present in the prioritized candidate variant table, using the Ensembl
REST API. One FASTA file per protein is written to --output, named
<GENE>_<ENSP>.fasta, ready for AlphaFold/ColabFold/SWISS-MODEL input
and for mutate_sequence.py.

Requires internet access to https://rest.ensembl.org (or set
ENSEMBL_REST_URL env var to a mirror/local instance).
"""

import argparse
import csv
import os
import sys
import time
import urllib.request
import urllib.error

ENSEMBL_REST_URL = os.environ.get("ENSEMBL_REST_URL", "https://rest.ensembl.org")


def fetch_protein_fasta(ensp_id, retries=3):
    url = f"{ENSEMBL_REST_URL}/sequence/id/{ensp_id}?type=protein;content-type=text/x-fasta"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(2 * (attempt + 1))
                continue
            sys.stderr.write(f"WARNING: HTTP {e.code} fetching {ensp_id}: {e}\n")
            return None
        except urllib.error.URLError as e:
            sys.stderr.write(f"WARNING: could not reach Ensembl REST for {ensp_id}: {e}\n")
            return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True, help="TSV from prioritize_variants.py")
    ap.add_argument("--output", required=True, help="output directory for wild-type FASTA files")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    seen = set()
    fetched = 0
    skipped = 0

    with open(args.variants) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ensp = row.get("ENSP", "").strip()
            gene = row.get("GENE", "unknown").strip()

            if not ensp or ensp == ".":
                sys.stderr.write(f"WARNING: no ENSP for {gene} {row.get('HGVSp','')} — skipping\n")
                skipped += 1
                continue

            key = (gene, ensp)
            if key in seen:
                continue
            seen.add(key)

            out_path = os.path.join(args.output, f"{gene}_{ensp}.fasta")
            if os.path.exists(out_path):
                continue

            fasta = fetch_protein_fasta(ensp)
            if fasta is None:
                skipped += 1
                continue

            with open(out_path, "w") as out:
                out.write(fasta)
            fetched += 1
            time.sleep(0.34)  # respect Ensembl REST rate limit (~3 req/s)

    sys.stderr.write(f"Fetched {fetched} wild-type protein sequences, skipped {skipped}, into {args.output}\n")


if __name__ == "__main__":
    main()
