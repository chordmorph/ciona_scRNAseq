#!/bin/bash

working_dir=
read1=
read2=
sample_name=
threads_num=
num_mismatches=
num_barcodes=
read_group=
library_name=

cd $working_dir
mkdir ./tmp/
mkdir ./tmp/sequence_files
mkdir ./tmp/summary_outputs
mkdir ./split_bams
mkdir ./split_bams/resampled_bams

#make reference files
#make a fasta index file for the genome fasta

samtools faidx ./genome_files/genome.fa

STAR \
  --runThreadN 32 \
  --runMode genomeGenerate \
  --genomeDir ./genome_files \
  --genomeFastaFiles ./genome_files/genome.fa \
  --sjdbGTFfile ./genome_files/annotation.gtf \
  --sjdbOverhang 115

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=genome_files/genome.fa \
  O=genome_files/genome.dict

gff2bed < ./genome_files/SNPS.gff3 > ./genome_files/ciona_snps.bed
gatk IndexFeatureFile -I ./genome_files/ciona_snps.bed

#convert the fastq reads to a paired read bam file
java -jar $EBROOTPICARD/picard.jar FastqToSam \
  F1=$read1 \
  F2=$read2 \
  OUTPUT=tmp/sequence_files/unaligned.bam \
  SM=$sample_name \
  PL=ILLUMINA \
  RG=$read_group \
  LB=$library_name

#tag the genome read with the barcodes and filter the reads (from dropseqtools)

TagBamWithReadSequenceExtended \
  INPUT=tmp/sequence_files/unaligned.bam \
  OUTPUT=tmp/sequence_files/unaligned_tagged_Cell.bam \
  SUMMARY=/tmp/summary_outputs/unaligned_tagged_Cellular.bam_summary.txt \
  BASE_RANGE=1-12 \
  BASE_QUALITY=10 \
  BARCODED_READ=1 \
  DISCARD_READ=False \
  TAG_NAME=XC \
  NUM_BASES_BELOW_QUALITY=1

TagBamWithReadSequenceExtended \
  INPUT=tmp/sequence_files/unaligned_tagged_Cell.bam \
  OUTPUT=tmp/sequence_files/unaligned_tagged.bam \
  SUMMARY=/tmp/summary_outputs/unaligned_tagged_Molecular.bam_summary.txt \
  BASE_RANGE=13-20 \
  BASE_QUALITY=10 \
  BARCODED_READ=1 \
  DISCARD_READ=True \
  TAG_NAME=XM \
  NUM_BASES_BELOW_QUALITY=1

FilterBam TAG_REJECT=XQ \
  INPUT=tmp/sequence_files/unaligned_tagged.bam \
  OUTPUT=tmp/sequence_files/unaligned_tagged_filtered.bam

#convert the bam file to a fastq for star Alignment
java -jar $EBROOTPICARD/picard.jar SamToFastq \
  INPUT=tmp/sequence_files/unaligned_tagged_filtered.bam \
  FASTQ=tmp/sequence_files/unaligned_tagged_filtered.fastq


#align the reads using a two pass alignment with STAR
STAR \
  --runThreadN $threads_num \
  --genomeDir ./genome_files \
  --readFilesIn tmp/sequence_files/unaligned_tagged_filtered.fastq \
  --outFilterMismatchNmax $num_mismatches \
  --twopassMode Basic \
  --outFileNamePrefix tmp/sequence_files/

#sort the aligned sam file by query order and save as a bam
java -jar $EBROOTPICARD/picard.jar SortSam \
  I=tmp/sequence_files/Aligned.out.sam \
  O=tmp/sequence_files/aligned.sorted.bam \
  SO=queryname

#merge the bam alignment to recover the cell and molecular tags

java -jar $EBROOTPICARD/picard.jar MergeBamAlignment \
  ALIGNED=tmp/sequence_files/aligned.sorted.bam \
  UNMAPPED=tmp/sequence_files/unaligned_tagged_filtered.bam \
  O=tmp/sequence_files/aligned.merged.bam \
  R=genome_files/genome.fa

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I=tmp/sequence_files/aligned.merged.bam \
  O=tmp/sequence_files/marked_duplicates.bam \
  M=tmp/summary_outputs/marked_dup_metrics.txt \
  BARCODE_TAG=XC

#clean up bead errors to match barcodes with dropseqpipe output

DetectBeadSynthesisErrors \
  I=tmp/sequence_files/marked_duplicates.bam \
  O=tmp/sequence_files/clean_marked_duplicates.bam \
  REPORT=tmp/summary_outputs/my_clean.indel_report.txt \
  OUTPUT_STATS=tmp/summary_outputs/my.synthesis_stats.txt \
  SUMMARY=tmp/summary_outputs/my.synthesis_stats.summary.txt \
  PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
  NUM_THREADS=$threads_num

gatk SplitNCigarReads \
  -R ./genome_files/genome.fa \
  -I tmp/sequence_files/clean_marked_duplicates.bam \
  -O tmp/sequence_files/merged_aligned_split.bam

gatk BaseRecalibrator \
  -I tmp/sequence_files/merged_aligned_split.bam \
  -R genome_files/genome.fa \
  -O tmp/sequence_files/recal_data.table \
  --known-sites genome_files/ciona_snps.bed

gatk ApplyBQSR \
  -R genome_files/genome.fa \
  -I tmp/sequence_files/merged_aligned_split.bam \
  --bqsr-recal-file tmp/sequence_files/recal_data.table \
  -O tmp/sequence_files/recalibrated_reads.bam

gatk AnalyzeCovariates \
  -bqsr tmp/sequence_files/recal_data.table \
  -plots tmp/summary_outputs/AnalyzeCovariates.pdf

#split bam by barcodes for variant calling
bamCleave \
  -b tmp/sequence_files/recalibrated_reads.bam \
  -o split_bams/$sample_name \
  -t XC \
  -c $num_barcodes

for bam_file in split_bams/*_sel_*.bam
do
  barcode=$(echo ${bam_file#*_sel_*} | cut -d '.' -f 1)
  full_sample_name=$(printf "%s_%s" $sample_name $barcode)
  samtools view -H $bam_file | sed "s/$sample_name/$full_sample_name/" | samtools reheader - $bam_file > split_bams/resampled_bams/$full_sample_name.bam
done
