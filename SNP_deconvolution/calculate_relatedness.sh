#!/bin/bash

SAMPLE=
SAMPLE_PREFIX=

vcftools \
  --vcf ${SAMPLE_PREFIX}/outputs/${SAMPLE}_raw_variants.vcf \
  --relatedness \
  --out ${SAMPLE_PREFIX}/outputs/${SAMPLE}_raw_variants_all_sites

vcftools \
  --vcf ${SAMPLE_PREFIX}/outputs/${SAMPLE}_raw_variants.vcf \
  --min-meanDP 0.5 \
  --relatedness \
  --out ${SAMPLE_PREFIX}/outputs/${SAMPLE}_raw_variants_subset_sites
