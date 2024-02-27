#!/bin/bash
# Author: Saul Pierotti, Birney reasearch, European Bioinformatics Institute

# Fills missing genotype calls in a VCF file from the DS field using a back and forth conversion with plink2
# Useful for the output from the birneylab/stitchimpute pipeline

# Edit parameters as required
N_CHR=24 # 24 for medaka
VCF="$1" # vcf from the command line
DOSAGE_FIELD="DS"
MEM="10000"
N_CPU="10"

plink2 \
  --memory $MEM \
  --threads $N_CPU \
  --out tmp \
  --chr-set $N_CHR \
  --vcf $VCF dosage=$DOSAGE_FIELD \
  --make-pgen vzs fill-missing-from-dosage erase-dosage

plink2 \
  --memory $MEM \
  --threads $N_CPU \
  --out "$(basename $VCF).no_missing"\
  --chr-set $N_CHR \
  --pgen tmp.pgen \
  --psam tmp.psam \
  --pvar tmp.var.zst \
  --export vcf bgz
