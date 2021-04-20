#!/bin/bash

working_dir=
sample1_split_dir=
sample2_split_dir=

ulimit -Sn 15000

cd $working_dir

mkdir ./tmp/
mkdir ./tmp/sequence_files/
mkdir ./tmp/split_raw_vcfs/
mkdir ./tmp/split_filtered_vcfs/
mkdir ./tmp/varclust_input_vcfs/
mkdir ./tmp/varclust_profiles/
mkdir ./outputs

samtools faidx ./genome_files/genome.fa

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=genome_files/genome.fa \
  O=genome_files/genome.dict

echo $(printf -- '-I %s ' $sample1_split_dir/*.bam) > tmp/input_file1.txt
echo $(printf -- '-I %s ' $sample2_split_dir/*.bam) > tmp/input_file2.txt
echo -O tmp/sequence_files/merged_multisample.bam > tmp/output_file.txt

gatk MergeSamFiles \
      --arguments_file tmp/input_file1.txt \
      --arguments_file tmp/input_file2.txt \
      --arguments_file tmp/output_file.txt

samtools index tmp/sequence_files/merged_multisample.bam
