#!/usr/bin/env python3
"""
mutate_sequence.py

Reads the prioritized candidate variant table plus the wild-type FASTA
files fetched by extract_protein_sequences.py, and generates mutant
protein FASTA files for missense variants by substituting the residue
named in HGVSp.

Supports HGVSp formats such as:
    p.Arg217Cys
    p.(Arg217Cys)
    ENSP00000384371.2:p.(Ala493Thr)

Non-missense HGVSp notations (frameshift, stop-gain, splice, etc.) are
logged and skipped, since a simple single-residue substitution does not
apply.

IMPORTANT:
The script verifies the reference residue named in HGVSp against the
fetched wild-type sequence at that position before writing a mutant.
If the residue does not match, the mutant is NOT written.

This HGVSp/sequence consistency check must pass before trusting the
variant for downstream FoldX analysis.
"""

import argparse
import csv
import os
import re
import sys


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
    "Ter": "*",
}


# Accept:
#   p.Ala493Thr
#   p.(Ala493Thr)
#   p.Arg217Cys
#   p.(Arg217Cys)
#
# The parentheses are optional because VEP commonly outputs:
#   p.(Ala493Thr)
HGVSP_MISSENSE_RE = re.compile(
    r"^p\.\(?(?P<ref>[A-Za-z]{3})(?P<pos>\d+)(?P<alt>[A-Za-z]{3})\)?$"
)


def parse_fasta(path):
    """Read a FASTA file and return (header, sequence)."""

    header = None
    seq = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                header = line
            else:
                seq.append(line)

    return header, "".join(seq)


def normalize_hgvsp(hgvsp):
    """
    Extract the protein-level HGVSp notation from a VEP field.

    Examples:
        ENSP00000384371.2:p.(Ala493Thr)
        ENST00000384371.7:c.1477G>A
        p.(Ala493Thr)

    Returns the final protein notation when available.
    """

    if not hgvsp:
        return ""

    # VEP fields may contain transcript/ENSP prefixes separated by ':'
    notation = hgvsp.split(":")[-1].strip()

    return notation


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--variants",
        required=True,
        help="Prioritized candidate variant TSV"
    )

    ap.add_argument(
        "--wildtype",
        required=True,
        help="Directory containing WT FASTA files"
    )

    ap.add_argument(
        "--output",
        required=True,
        help="Output directory for mutant FASTA files"
    )

    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    n_written = 0
    n_skipped = 0
    n_mismatch = 0

    with open(args.variants) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            gene = row.get("GENE", "").strip()
            ensp = row.get("ENSP", "").strip()
            hgvsp = row.get("HGVSp", "").strip()

            # --------------------------------------------------------
            # Basic validation
            # --------------------------------------------------------

            if (
                not gene
                or gene == "."
                or not ensp
                or ensp == "."
                or not hgvsp
                or hgvsp == "."
            ):
                n_skipped += 1
                continue

            # --------------------------------------------------------
            # Normalize HGVSp
            # --------------------------------------------------------

            hgvsp_notation = normalize_hgvsp(hgvsp)

            # --------------------------------------------------------
            # Parse missense variant
            # --------------------------------------------------------

            match = HGVSP_MISSENSE_RE.match(hgvsp_notation)

            if not match:

                sys.stderr.write(
                    f"SKIP (not a simple missense HGVSp): "
                    f"{gene} {hgvsp_notation}\n"
                )

                n_skipped += 1
                continue

            ref_aa3 = match.group("ref")
            pos = int(match.group("pos"))
            alt_aa3 = match.group("alt")

            ref_aa1 = AA3_TO_1.get(ref_aa3)
            alt_aa1 = AA3_TO_1.get(alt_aa3)

            # --------------------------------------------------------
            # Validate amino-acid codes
            # --------------------------------------------------------

            if (
                ref_aa1 is None
                or alt_aa1 is None
                or alt_aa1 == "*"
            ):

                sys.stderr.write(
                    f"SKIP (unsupported residue code): "
                    f"{gene} {hgvsp_notation}\n"
                )

                n_skipped += 1
                continue

            # --------------------------------------------------------
            # Locate WT FASTA
            # --------------------------------------------------------

            wt_path = os.path.join(
                args.wildtype,
                f"{gene}_{ensp}.fasta"
            )

            if not os.path.exists(wt_path):

                sys.stderr.write(
                    f"SKIP (no WT fasta found): "
                    f"{wt_path}\n"
                )

                n_skipped += 1
                continue

            # --------------------------------------------------------
            # Read WT sequence
            # --------------------------------------------------------

            header, seq = parse_fasta(wt_path)

            if not seq:

                sys.stderr.write(
                    f"SKIP (empty WT sequence): "
                    f"{wt_path}\n"
                )

                n_skipped += 1
                continue

            # --------------------------------------------------------
            # Check position is within sequence
            # --------------------------------------------------------

            if pos < 1 or pos > len(seq):

                sys.stderr.write(
                    f"MISMATCH: {gene} {hgvsp_notation} "
                    f"position {pos} outside sequence length "
                    f"{len(seq)} — check HGVSp/transcript build match\n"
                )

                n_mismatch += 1
                continue

            # --------------------------------------------------------
            # Verify reference amino acid
            # --------------------------------------------------------

            observed_aa = seq[pos - 1]

            if observed_aa != ref_aa1:

                sys.stderr.write(
                    f"MISMATCH: {gene} {hgvsp_notation} "
                    f"expected {ref_aa1} at position {pos} "
                    f"but WT sequence has {observed_aa} — "
                    f"verify HGVSp field before trusting this variant. "
                    f"Mutant NOT written.\n"
                )

                n_mismatch += 1
                continue

            # --------------------------------------------------------
            # Generate mutant sequence
            # --------------------------------------------------------

            mutant_seq = (
                seq[:pos - 1]
                + alt_aa1
                + seq[pos:]
            )

            mutant_id = (
                f"{gene}_{ensp}_{ref_aa1}{pos}{alt_aa1}"
            )

            out_path = os.path.join(
                args.output,
                f"{mutant_id}.fasta"
            )

            # --------------------------------------------------------
            # Write mutant FASTA
            # --------------------------------------------------------

            with open(out_path, "w") as out:

                out.write(
                    f">{mutant_id} | "
                    f"WT:{header.lstrip('>')} | "
                    f"variant:{hgvsp_notation}\n"
                )

                for i in range(0, len(mutant_seq), 60):

                    out.write(
                        mutant_seq[i:i + 60] + "\n"
                    )

            n_written += 1

    # ------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------

    sys.stderr.write(
        f"Mutant sequences written: {n_written}, "
        f"skipped: {n_skipped}, "
        f"HGVSp/sequence mismatches: {n_mismatch}\n"
    )

    if n_mismatch:

        sys.stderr.write(
            "\n"
            "WARNING: HGVSp/sequence mismatches were detected.\n"
            "These usually indicate that the ENSP sequence used by "
            "VEP differs from the WT sequence being used here, or "
            "that the HGVSp annotation needs verification.\n"
            "Mismatching variants were NOT written and should not "
            "be used for FoldX until resolved.\n"
        )


if __name__ == "__main__":
    main()

