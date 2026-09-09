#!/usr/bin/env python3

"""
Convert prioritized HGVS protein variants into FoldX individual_list.txt format.

Example
-------
Input:
    ENSP00000269305:p.Arg175His

Output:
    RA175H;

FoldX format:
WT_CHAIN_POSITION_MUT;
"""

import argparse
import re
import sys

import pandas as pd
from Bio import SeqIO


AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}

HGVS3 = re.compile(
    r".*p\.\(?(?P<wt>[A-Za-z]{3})(?P<pos>\d+)(?P<mut>[A-Za-z]{3})\)?"
)

HGVS1 = re.compile(
    r".*p\.\(?(?P<wt>[A-Za-z])(?P<pos>\d+)(?P<mut>[A-Za-z])\)?"
)


def aa_to_one(code):

    code = code.strip()

    if len(code) == 1:
        return code.upper()

    return AA3_TO_1.get(code.capitalize())


def parse_hgvs(hgvs):

    hgvs = str(hgvs)

    if any(x in hgvs for x in ["fs", "*", "="]):
        return None

    m = HGVS3.match(hgvs)

    if m:

        wt = aa_to_one(m.group("wt"))
        mut = aa_to_one(m.group("mut"))

        if wt and mut:
            return wt, m.group("pos"), mut

    m = HGVS1.match(hgvs)

    if m:

        return (
            m.group("wt").upper(),
            m.group("pos"),
            m.group("mut").upper(),
        )

    return None


def load_sequence(fasta):

    seq = ""

    for record in SeqIO.parse(fasta, "fasta"):
        seq += str(record.seq)

    return seq


def validate(sequence, pos, wt):

    try:
        return sequence[int(pos) - 1] == wt
    except IndexError:
        return False


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--variants", required=True)
    parser.add_argument("--canonical-map", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--gene-col",
        default="GENE"
    )

    parser.add_argument(
        "--protein-col",
        default="HGVSp"
    )

    parser.add_argument(
        "--chain",
        default="A"
    )

    parser.add_argument(
        "--fasta",
        default=None
    )

    args = parser.parse_args()

    ############################################################
    # Load files
    ############################################################

    df = pd.read_csv(
        args.variants,
        sep="\t"
    )

    canonical = (
        pd.read_csv(
            args.canonical_map,
            sep="\t"
        )
        .set_index("gene")["ensp"]
        .to_dict()
    )

    ############################################################
    # Validate columns
    ############################################################

    required = [
        args.gene_col,
        args.protein_col,
        "ENSP",
        "CONSEQUENCE"
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        sys.exit(
            f"Missing required columns: {', '.join(missing)}"
        )

    ############################################################
    # Check canonical transcript exists
    ############################################################

    if args.gene not in canonical:

        sys.exit(
            f"{args.gene} not found in canonical ENSP mapping."
        )

    ############################################################
    # Keep only canonical missense variants
    ############################################################

    gene_df = df[
        (df[args.gene_col] == args.gene)
        &
        (df["ENSP"] == canonical[args.gene])
        &
        (
            df["CONSEQUENCE"]
            .fillna("")
            .str.contains("missense_variant")
        )
    ]

    ############################################################
    # No compatible variants
    ############################################################

    if gene_df.empty:

        print(
            f"No FoldX-compatible variants for {args.gene}"
        )

        open(args.output, "w").close()

        return

    ############################################################
    # Optional sequence validation
    ############################################################

    sequence = None

    if args.fasta:

        sequence = load_sequence(args.fasta)

        print(
            f"Loaded sequence ({len(sequence)} aa)"
        )

    ############################################################
    # Convert HGVS
    ############################################################

    mutations = []

    skipped = []

    mismatches = []

    for hgvs in gene_df[
        args.protein_col
    ].dropna().unique():

        parsed = parse_hgvs(hgvs)

        if parsed is None:

            skipped.append(hgvs)

            continue

        wt, pos, mut = parsed

        if sequence:

            if not validate(
                sequence,
                pos,
                wt
            ):

                mismatches.append(
                    hgvs
                )

                continue

        mutations.append(
            f"{wt}{args.chain}{pos}{mut};"
        )

    ############################################################
    # Remove duplicates
    ############################################################

    mutations = sorted(
        set(mutations),
        key=lambda x: int(
            re.search(r"\d+", x).group()
        )
    )

    ############################################################
    # Logging
    ############################################################

    if skipped:

        print(
            "\nSkipped HGVS:",
            file=sys.stderr
        )

        for x in skipped:
            print(
                x,
                file=sys.stderr
            )

    if mismatches:

        print(
            "\nSequence mismatches:",
            file=sys.stderr
        )

        for x in mismatches:
            print(
                x,
                file=sys.stderr
            )

    ############################################################
    # No valid mutations
    ############################################################

    if not mutations:

        print(
            f"No valid FoldX mutations for {args.gene}"
        )

        open(args.output, "w").close()

        return

    ############################################################
    # Write output
    ############################################################

    with open(
        args.output,
        "w"
    ) as out:

        out.write(
            "\n".join(mutations)
        )

        out.write("\n")

    print(
        f"Generated {len(mutations)} mutation(s)"
    )

    print(
        f"Output: {args.output}"
    )


if __name__ == "__main__":
    main()
