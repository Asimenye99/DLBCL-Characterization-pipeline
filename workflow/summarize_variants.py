import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep="\t", header=None)

df.columns = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Gene",
    "Consequence",
    "Impact",
    "HGVSc",
    "HGVSp",
    "SIFT",
    "PolyPhen"
]

print("\nTotal variants:", len(df))

print("\nImpact counts")
print(df["Impact"].value_counts())

print("\nTop mutated genes")
print(df["Gene"].value_counts().head(20))
