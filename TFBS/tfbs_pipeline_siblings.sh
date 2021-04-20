#!/bin/bash

MARKERS_PATH=#directory with text files of target and control gene names
CLUSTER_NAMES=#path to text file with names of clusters to analyze, must be identical to names used in markers directory

FASTA_OUTPUT_DIR=#directory path for target and control sequence input to opossum
OPOSSUM_OUTPUT_DIR=#directory path for opossum output

PWM_PATH=#path to JASPAR pwm file to use
FASTA_PATH=#path to genome fasta
FASTA_INDEX_PATH=#path to genome fasta index
GFF_PATH=#path to gene GFF
ATAC_PEAKS_BED_PATH=#path to ATAC peaks BED file

for CLUSTER_NAME in $(cat ${CLUSTER_NAMES})
do
  [  -d ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME} ] && rm -rf ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}
  [  -d ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME} ] && rm -rf ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME}
  mkdir ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}
  mkdir ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME}

  submit 1 \
    -N MAKE_${CLUSTER_NAME}_fastas \
    python tfbs_enrichment_analysis.py \
      ${MARKERS_PATH}/${CLUSTER_NAME}_target_gene_models.csv \
      ${MARKERS_PATH}/${CLUSTER_NAME}_control_gene_models.csv \
      ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}/ \
      1500 \
      up \
      False \
      ${GFF_PATH} \
      ${FASTA_INDEX_PATH} \
      ${FASTA_PATH} \
      ${ATAC_PEAKS_BED_PATH} \
      yes

  submit 1 \
    -N RUN_OPOSSUM_${CLUSTER_NAME}_NO_ATAC \
    -hold_jid MAKE_${CLUSTER_NAME}_fastas \
    opossum_seq_ssa.pl \
      -d ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME}/opossum_output \
      -s ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}/target_sequences.fa \
      -b ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}/control_sequences.fa \
      -mf ${PWM_PATH} \
      -nh \
      -sr zscore

  submit 1 \
    -N RUN_OPOSSUM_${CLUSTER_NAME}_WITH_ATAC \
    -hold_jid MAKE_${CLUSTER_NAME}_fastas \
    opossum_seq_ssa.pl \
      -d ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME}/opossum_output_ATAC \
      -s ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}/target_sequences_ATAC_peaks.fa \
      -b ${FASTA_OUTPUT_DIR}/${CLUSTER_NAME}/control_sequences_ATAC_peaks.fa \
      -mf ${PWM_PATH} \
      -nh \
      -sr zscore

  submit 1 \
    -N AVERAGE_${CLUSTER_NAME}_RESULTS \
    -hold_jid RUN_OPOSSUM_${CLUSTER_NAME}_WITH_ATAC,RUN_OPOSSUM_${CLUSTER_NAME}_NO_ATAC \
    python averaging_script.py \
      ${OPOSSUM_OUTPUT_DIR}/${CLUSTER_NAME} \
      ${CLUSTER_NAME}
done
