#!/usr/bin/env python3
"""
mutate_sequence.py

Reads the prioritized candidate variant table plus the wild-type FASTA
files fetched by extract_protein_sequences.py, and generates mutant
protein FASTA files (.fasta, matching the wild-type naming convention)
for missense variants by substituting the residue named in HGVSp
(e.g. p.Arg217Cys / p.R217C).

Non-missense HGVSp notations (frameshift, stop-gain, splice, etc.) are
logged and skipped, since a simple single-residue substitution does not
apply.

IMPORTANT: this script verifies the reference residue named in HGVSp
against the fetched wild-type sequence at that position before writing
a mutant, and will refuse — with a clear warning — to write anything if
there is a mismatch. This is the HGVSp/sequence consistency check that
needs to pass before FoldX input.
"""

import argparse
import csv
import os
import re
import sys

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

HGVSP_MISSENSE_RE = re.compile(
    r"p\.(?P<ref>[A-Za-z]{3})(?P<pos>\d+)(?P<alt>[A-Za-z]{3})$"
)


def parse_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
            else:
                seq.append(line)
    return header, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--wildtype", required=True, help="directory of WT FASTA files from extract_protein_sequences.py")
    ap.add_argument("--output", required=True, help="output directory for mutant FASTA files")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    n_written = 0
    n_skipped = 0
    n_mismatch = 0

    with open(args.variants) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("GENE", "unknown").strip()
            ensp = row.get("ENSP", "").strip()
            hgvsp = row.get("HGVSp", "").strip()

            if not ensp or ensp == "." or not hgvsp or hgvsp == ".":
                n_skipped += 1
                continue

            hgvsp_notation = hgvsp.split(":")[-1]
            m = HGVSP_MISSENSE_RE.match(hgvsp_notation)
            if not m:
                sys.stderr.write(
                    f"SKIP (not a simple missense HGVSp): {gene} {hgvsp_notation}\n"
                )
                n_skipped += 1
                continue

            ref_aa3, pos, alt_aa3 = m.group("ref"), int(m.group("pos")), m.group("alt")
            ref_aa1 = AA3_TO_1.get(ref_aa3)
            alt_aa1 = AA3_TO_1.get(alt_aa3)

            if ref_aa1 is None or alt_aa1 is None or alt_aa1 == "*":
                sys.stderr.write(f"SKIP (unsupported residue code): {gene} {hgvsp_notation}\n")
                n_skipped += 1
                continue

            wt_path = os.path.join(args.wildtype, f"{gene}_{ensp}.fasta")
            if not os.path.exists(wt_path):
                sys.stderr.write(f"SKIP (no WT fasta found): {wt_path}\n")
                n_skipped += 1
                continue

            header, seq = parse_fasta(wt_path)

            if pos < 1 or pos > len(seq):
                sys.stderr.write(
                    f"MISMATCH: {gene} {hgvsp_notation} position {pos} outside "
                    f"sequence length {len(seq)} — check HGVSp/transcript build match\n"
                )
                n_mismatch += 1
                continue

            observed_aa = seq[pos - 1]
            if observed_aa != ref_aa1:
                sys.stderr.write(
                    f"MISMATCH: {gene} {hgvsp_notation} expected {ref_aa1} at position "
                    f"{pos} but WT sequence has {observed_aa} — verify HGVSp field "
                    f"before trusting this variant. Mutant NOT written.\n"
                )
                n_mismatch += 1
                continue

            mutant_seq = seq[:pos - 1] + alt_aa1 + seq[pos:]
            mutant_id = f"{gene}_{ensp}_{ref_aa1}{pos}{alt_aa1}"
            out_path = os.path.join(args.output, f"{mutant_id}.fasta")

            with open(out_path, "w") as out:
                out.write(f">{mutant_id} | WT:{header.lstrip('>')} | variant:{hgvsp_notation}\n")
                for i in range(0, len(mutant_seq), 60):
                    out.write(mutant_seq[i:i + 60] + "\n")

            n_written += 1

    sys.stderr.write(
        f"Mutant sequences written: {n_written}, skipped: {n_skipped}, "
        f"HGVSp/sequence mismatches: {n_mismatch}\n"
    )
    if n_mismatch:
        sys.stderr.write(
            "NOTE: mismatches usually mean the ENSP build used by VEP's cache "
            "differs from the REST API's default sequence, or HGVSp needs "
            "re-verification, before trusting this variant for FoldX input.\n"
        )


if __name__ == "__main__":
    main()
