#!/usr/bin/env python3
"""Remove FoldX mutations whose wild-type residue is absent from the structure.

A mutation such as PA414L is kept only if chain A of the structure has
proline at position 414. FoldX stops with 'Specified residue not found' when
this is not the case, and a single such entry makes the whole gene fail. The
usual cause is that VEP numbers the variant on the Ensembl canonical protein
while AlphaFold uses the UniProt canonical sequence (different isoform, or a
position beyond the end of the UniProt protein).

Removed mutations are written to --unmodelled with the reason, so they are
reported as 'not modelled' rather than lost.
"""
import argparse, os

T = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
     "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
     "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", required=True)
ap.add_argument("--mutlist", required=True, help="FoldX individual_list (edited in place)")
ap.add_argument("--unmodelled", required=True)
ap.add_argument("--gene", default="")
a = ap.parse_args()

res = {}   # (chain, resnum) -> one-letter
if os.path.exists(a.pdb) and os.path.getsize(a.pdb) > 0:
    for l in open(a.pdb):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            res[(l[21], l[22:27].strip())] = T.get(l[17:20].strip(), "X")

keep, drop = [], []
for line in open(a.mutlist):
    line = line.strip()
    if not line:
        continue
    reasons = []
    for tok in [t for t in line.split(";") if t]:
        wt, ch, pos = tok[0], tok[1], tok[2:-1]
        have = res.get((ch, pos))
        if not res:
            reasons.append("no AlphaFold structure")
        elif have is None:
            reasons.append(f"position {pos} outside structure (protein length {max(int(p) for _, p in res if p.isdigit())})")
        elif have != wt:
            reasons.append(f"residue {pos} is {have} in structure, {wt} expected")
    (drop if reasons else keep).append((line, "; ".join(reasons)))

with open(a.mutlist, "w") as f:
    for line, _ in keep:
        f.write(line + "\n")
with open(a.unmodelled, "w") as f:
    f.write("gene\tmutation\treason\n")
    for line, why in drop:
        f.write(f"{a.gene}\t{line}\t{why}\n")
print(f"{a.gene}: kept {len(keep)}, not modelled {len(drop)}")
