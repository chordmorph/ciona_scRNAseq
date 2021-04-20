#!/bin/bash

SAMPLE_PREFIX=

gatk HaplotypeCaller  \
  -R genome_files/genome.fa \
  -I ${SAMPLE_PREFIX}/tmp/sequence_files/merged_multisample.bam \
  -O ${SAMPLE_PREFIX}/tmp/sequence_files/raw_variants_multisample.vcf \
  --disable-read-filter NotDuplicateReadFilter
