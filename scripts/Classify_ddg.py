#!/usr/bin/env python3

"""
Parse FoldX BuildModel output (Average_*.fxout) for each gene
and classify the structural stability effect of each mutation.

For each gene:
- Reads FoldX Average_*.fxout
- Reads Gene_individual_list.txt
- Extracts FoldX averaged ΔΔG values
- Classifies variants as:
    * destabilizing (>0.5 kcal/mol)
    * stabilizing (<-0.5 kcal/mol)
    * neutral (-0.5 to 0.5 kcal/mol)

Outputs a TSV summarizing the structural effect of each variant.
"""

import argparse
import glob
import os
import sys

import pandas as pd


def read_avg_fxout(path):
    """
    Read FoldX Average_*.fxout.

    Returns a dataframe with:
        Pdb
        ddG
    """

    rows = []

    with open(path) as fh:
        for line in fh:

            line = line.strip()

            if (
                not line
                or line.startswith("FoldX")
                or line.startswith("by")
                or line.startswith("Jesper")
                or line.startswith("Luis")
                or line.startswith("-")
                or line.startswith("PDB file")
                or line.startswith("Output type")
                or line.startswith("Pdb")
            ):
                continue

            fields = line.split()

            if len(fields) < 3:
                continue

            try:
                ddg = float(fields[2])
            except ValueError:
                continue

            rows.append(
                {
                    "Pdb": fields[0],
                    "ddG": ddg,
                }
            )

    return pd.DataFrame(rows)


def read_mutation_list(path):
    """
    Read FoldX Gene_individual_list.txt.

    Removes the trailing semicolon.
    """

    with open(path) as fh:
        return [
            line.strip().rstrip(";")
            for line in fh
            if line.strip()
        ]


def classify(ddg, destabilizing_threshold, stabilizing_threshold):
    """
    Classify stability change from ΔΔG.
    """

    if ddg > destabilizing_threshold:
        return "destabilizing"

    if ddg < stabilizing_threshold:
        return "stabilizing"

    return "neutral"


def main():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--buildmodel-dir",
        required=True,
        help="Directory containing FoldX BuildModel output, e.g. results/foldx/buildmodel/"
    )

    parser.add_argument(
        "--mutation-dir",
        required=True,
        help="Directory containing Gene_individual_list.txt files, e.g. results/foldx/mutation_lists/"
    )

    parser.add_argument(
        "--genes",
        nargs="+",
        required=True
    )

    parser.add_argument(
        "--destabilizing-threshold",
        type=float,
        default=0.5
    )

    parser.add_argument(
        "--stabilizing-threshold",
        type=float,
        default=-0.5
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    all_rows = []

    for gene in args.genes:

        print(f"Processing {gene}...", file=sys.stderr)

        # ---------------------------------------------------------
        # FoldX BuildModel directory
        # ---------------------------------------------------------

        gene_dir = os.path.join(
            args.buildmodel_dir,
            gene
        )

        # ---------------------------------------------------------
        # Mutation list
        #
        # Example:
        # results/foldx/mutation_lists/MEF2B_individual_list.txt
        # ---------------------------------------------------------

        mutlist_path = os.path.join(
            args.mutation_dir,
            f"{gene}_individual_list.txt"
        )

        # ---------------------------------------------------------
        # Locate Average FoldX output
        # ---------------------------------------------------------

        avg_candidates = sorted(
            glob.glob(
                os.path.join(
                    gene_dir,
                    "Average_*.fxout"
                )
            )
        )

        # ---------------------------------------------------------
        # Check mutation list
        # ---------------------------------------------------------

        if not os.path.exists(mutlist_path):

            print(
                f"WARNING: Missing mutation list for {gene}: "
                f"{mutlist_path}",
                file=sys.stderr
            )

            continue

        mutations = read_mutation_list(
            mutlist_path
        )

        if len(mutations) == 0:

            print(
                f"INFO: No FoldX-compatible mutations for {gene}. "
                f"Skipping.",
                file=sys.stderr
            )

            continue

        # ---------------------------------------------------------
        # Check FoldX output
        # ---------------------------------------------------------

        if not avg_candidates:

            print(
                f"WARNING: Missing FoldX Average output for {gene}: "
                f"{gene_dir}",
                file=sys.stderr
            )

            continue

        avg_path = avg_candidates[0]

        if os.path.getsize(avg_path) == 0:

            print(
                f"INFO: Empty FoldX output for {gene}. "
                f"Skipping.",
                file=sys.stderr
            )

            continue

        # ---------------------------------------------------------
        # Read FoldX output
        # ---------------------------------------------------------

        try:

            avg_df = read_avg_fxout(
                avg_path
            )

        except Exception as e:

            print(
                f"WARNING: Could not read FoldX output for "
                f"{gene}: {e}",
                file=sys.stderr
            )

            continue

        if avg_df.empty:

            print(
                f"WARNING: No parseable FoldX result rows for "
                f"{gene}",
                file=sys.stderr
            )

            continue

        avg_df = avg_df.reset_index(
            drop=True
        )

        # ---------------------------------------------------------
        # Check mutation/result counts
        # ---------------------------------------------------------

        if len(avg_df) != len(mutations):

            print(
                f"WARNING: {gene} has "
                f"{len(avg_df)} FoldX results but "
                f"{len(mutations)} mutations.",
                file=sys.stderr
            )

        # ---------------------------------------------------------
        # Match mutations to FoldX results
        # ---------------------------------------------------------

        for i, mutation in enumerate(mutations):

            if i >= len(avg_df):

                print(
                    f"WARNING: Missing FoldX result for "
                    f"{gene} {mutation}",
                    file=sys.stderr
                )

                continue

            ddg = avg_df.iloc[i]["ddG"]

            if pd.isna(ddg):
                continue

            all_rows.append(
                {
                    "gene": gene,
                    "variant": mutation,
                    "PDB": avg_df.iloc[i]["Pdb"],
                    "ddG_kcal_per_mol": round(
                        float(ddg),
                        3
                    ),
                    "classification": classify(
                        ddg,
                        args.destabilizing_threshold,
                        args.stabilizing_threshold,
                    ),
                }
            )

    # -------------------------------------------------------------
    # No results
    # -------------------------------------------------------------

    if not all_rows:

        print(
            "WARNING: No structural results were generated.",
            file=sys.stderr
        )

        pd.DataFrame(
            columns=[
                "gene",
                "variant",
                "PDB",
                "ddG_kcal_per_mol",
                "classification",
            ]
        ).to_csv(
            args.output,
            sep="\t",
            index=False,
        )

        print(
            f"Wrote empty results table to {args.output}"
        )

        return

    # -------------------------------------------------------------
    # Create output dataframe
    # -------------------------------------------------------------

    out_df = pd.DataFrame(
        all_rows,
        columns=[
            "gene",
            "variant",
            "PDB",
            "ddG_kcal_per_mol",
            "classification",
        ],
    )

    # Sort results
    out_df = out_df.sort_values(
        ["gene", "variant"]
    )

    # -------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------

    out_df.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    print(
        f"Wrote {len(out_df)} structural classifications "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
