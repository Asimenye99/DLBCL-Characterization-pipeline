#!/usr/bin/env python3

import argparse
import sys
import csv
import re

DLBCL_DRIVER_GENES = {
    "TP53", "MYD88", "CD79B", "CREBBP", "EZH2", "KMT2D", "BCL2", "BCL6",
    "TNFRSF14", "GNA13", "B2M", "CD58", "PIM1", "TBL1XR1", "SOCS1",
    "ITPKB", "ETS1", "MEF2B", "STAT6", "HIST1H1E", "TNFAIP3", "PRDM1",
    "DTX1", "ATM", "POU2F2", "IRF4", "BTG1", "BTG2", "FOXO1", "XPO1",
    "NOTCH2", "SGK1", "CCND3", "RHOA", "ZFP36L1", "S1PR2", "GNAI2",
    "P2RY8", "MYC", "PTEN", "TET2", "CARD11", "TCF3", "ID3", "BCL10",
    "NOTCH1", "NFKBIE",
}

HIGH_IMPACT = {"HIGH"}
MODERATE_IMPACT = {"MODERATE"}
CLINVAR_PATH_RE = re.compile(r"pathogenic", re.IGNORECASE)


def parse_amino_acids(aa, pos):
    """
    Primary FoldX source: Amino_acids field (e.g. R/G)
    """
    if not aa or "/" not in aa:
        return None, None, None, None

    if not pos:
        return None, None, None, None

    ref, alt = aa.split("/")
    ref = ref.strip()
    alt = alt.strip()
    short = f"p.{ref}{pos}{alt}"

    return ref, pos, alt, short


def parse_hgvsp(hgvsp):
    """
    Fallback HGVS protein parser
    """
    if not hgvsp:
        return None, None, None, None

    m = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*|=)", hgvsp)
    if not m:
        return None, None, None, None

    ref3, pos, alt3 = m.group(1), m.group(2), m.group(3)

    AA3TO1 = {
        "Ala": "A","Arg": "R","Asn": "N","Asp": "D","Cys": "C",
        "Gln": "Q","Glu": "E","Gly": "G","His": "H","Ile": "I",
        "Leu": "L","Lys": "K","Met": "M","Phe": "F","Pro": "P",
        "Ser": "S","Thr": "T","Trp": "W","Tyr": "Y","Val": "V",
    }

    ref1 = AA3TO1.get(ref3, ref3)
    alt1 = AA3TO1.get(alt3, alt3) if alt3 not in ["*", "="] else alt3

    short = f"p.{ref1}{pos}{alt1}"
    return ref1, pos, alt1, short


def main():
    ap = argparse.ArgumentParser(description="DLBCL driver gene protein extraction (FoldX-ready)")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="SAMPLE")
    args = ap.parse_args()

    from cyvcf2 import VCF

    vcf = VCF(args.vcf)

    # Parse CSQ header
    csq_header = None
    for h in vcf.header_iter():
        info = h.info(extra=True)
        if info.get("ID") == "CSQ":
            csq_header = info["Description"].split("Format:")[1].strip().strip('"').split("|")
            break

    if not csq_header:
        sys.exit("ERROR: CSQ header not found in VCF")

    idx = {k: i for i, k in enumerate(csq_header)}

    rows = []

    for variant in vcf:
        csq_raw = variant.INFO.get("CSQ")
        if not csq_raw:
            continue

        for ann in csq_raw.split(","):
            fields = ann.split("|")
            if len(fields) != len(csq_header):
                continue

            gene = fields[idx.get("SYMBOL", -1)]
            if gene.upper() not in DLBCL_DRIVER_GENES:
                continue

            consequence = fields[idx.get("Consequence", -1)]
            impact = fields[idx.get("IMPACT", -1)]
            clin_sig = fields[idx.get("CLIN_SIG", -1)] if "CLIN_SIG" in idx else ""
            transcript = fields[idx.get("Feature", -1)]

            # 🔥 STRICT FILTER: missense only
            if "missense_variant" not in consequence:
                continue

            # Optional but recommended: canonical only if present
            if "CANONICAL" in idx:
                canonical = fields[idx["CANONICAL"]]
                if canonical not in ["YES", "1"]:
                    continue

            aa = fields[idx.get("Amino_acids", -1)] if "Amino_acids" in idx else ""
            pos = fields[idx.get("Protein_position", -1)] if "Protein_position" in idx else ""
            hgvsp = fields[idx.get("HGVSp", -1)] if "HGVSp" in idx else ""

            # clinical filter (optional signal boost)
            is_path = bool(CLINVAR_PATH_RE.search(clin_sig)) if clin_sig else False

            if impact not in HIGH_IMPACT and impact not in MODERATE_IMPACT and not is_path:
                continue

            # PRIMARY: amino acids
            ref_aa, pos_aa, alt_aa, short = parse_amino_acids(aa, pos)

            # FALLBACK: HGVSp
            if not short:
                ref_aa, pos_aa, alt_aa, short = parse_hgvsp(hgvsp)

            if not short:
                continue

            rows.append({
                "SAMPLE": args.sample,
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": ",".join(variant.ALT),
                "GENE": gene,
                "TRANSCRIPT": transcript,
                "CONSEQUENCE": consequence,
                "IMPACT": impact,
                "CLIN_SIG": clin_sig,
                "HGVSp_raw": hgvsp,
                "AA_CHANGE": short,
                "REF_AA": ref_aa,
                "AA_POSITION": pos_aa,
                "ALT_AA": alt_aa,
            })

    rows.sort(key=lambda r: (r["GENE"], r["CHROM"], r["POS"]))

    fieldnames = [
        "SAMPLE","CHROM","POS","REF","ALT","GENE","TRANSCRIPT",
        "CONSEQUENCE","IMPACT","CLIN_SIG","HGVSp_raw",
        "AA_CHANGE","REF_AA","AA_POSITION","ALT_AA"
    ]

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    foldx = sum(1 for r in rows if r["AA_CHANGE"])

    print(f"Wrote {len(rows)} driver-gene variants to {args.out}")
    print(f"  -> {foldx} FoldX-ready AA changes")


if __name__ == "__main__":
    main()
