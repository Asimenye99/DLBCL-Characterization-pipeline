#!/bin/bash

INPUT=$1
OUTPUT=$2

bcftools +split-vep "$INPUT" \
-f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%HGVSc\t%HGVSp\t%SIFT\t%PolyPhen\n' \
> "$OUTPUT"
