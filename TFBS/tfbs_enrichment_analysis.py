#!/usr/bin/env python3

import pybedtools
import pandas as pd
import sys


def matcheslist(feature, list_to_match):
    gene = feature.attrs['ID']
    return gene in list_to_match

def matcheslist_mRNA(feature, list_to_match):
    gene = feature.attrs['Parent']
    return gene in list_to_match

def flank_filter_sort(list, bp_length, up_and_or_downstream, exclude_exons, exons, chrom_sizes):
    if up_and_or_downstream == 'up':
        list_flanked = list.flank(g=chrom_sizes,l=bp_length,r=0,s=True)
    else:
        list_flanked = list.flank(g=chrom_sizes,l=bp_length,r=0,s=True)
        list_flanked = list_flanked.slop(g=chrom_sizes, s=True, r=bp_length, l=0)
    if exclude_exons == True:
        list_flanked = list_flanked.subtract(exons).saveas()
    else:
        list_flanked = list_flanked
    filtered_flanked = list_flanked.filter(lambda x: (x.end-x.start)>100).saveas()
    sorted_filtered_flanked_list = filtered_flanked.sort().merge().saveas()
    return sorted_filtered_flanked_list



def tfbs_enrichment_analysis(genes_of_interest, genes_of_uninterest, output_dir,
                             bp_length, up_and_or_downstream, exclude_exons,
                             gene_model_path, chrom_sizes_path, genome_fa_path,
                             atac_bed_path, use_atac_peaks):

    #using GFF3 from aniseed (KH 2012 models with NCBI models (mito))
    ATACpeaks = pybedtools.BedTool(atac_bed_path)
    gene_models = pybedtools.BedTool(gene_model_path)
    allgenes = gene_models.filter(lambda x: x[2]=='gene').saveas()
    alloperons = gene_models.filter(lambda x: x[2]=='operon').saveas()
    alltranscripts = gene_models.filter(lambda x: x[2]=='mRNA').saveas()
    all_exons=gene_models.filter(lambda x: x[2]=='CDS').saveas()
    chromosmome_sizes = chrom_sizes_path
    genome = genome_fa_path

    #find enhancer regions of genes of interest
    up_genes_of_interest = pd.read_csv(genes_of_interest).gene_model.tolist()
    top_genes = allgenes.filter(matcheslist, list_to_match = up_genes_of_interest).saveas()

    #get regions of interest that lie in operons
    top_operons = alloperons.intersect(top_genes, u=True, s=True)

    #get regions of interest that are not in operons
    top_transcripts = alltranscripts.filter(matcheslist_mRNA, list_to_match=up_genes_of_interest).saveas()
    top_transcripts_no_operons = top_transcripts.subtract(alloperons, s=True, A=True)

    #find flanking regions for upregulated genes
    operons_upstream = top_operons.flank(g=chromosmome_sizes,l=bp_length,r=0,s=True)
    operons_upstream3 = operons_upstream.filter(lambda x: (x.end-x.start)>100).saveas()
    top_transcript_enhancers = flank_filter_sort(top_transcripts_no_operons, bp_length, up_and_or_downstream, exclude_exons, all_exons, chromosmome_sizes)
    top_enhancers = top_transcript_enhancers.cat(operons_upstream3, postmerge=True)



    #find enhancer regions of downregulated genes
    downregulated_genes = pd.read_csv(genes_of_uninterest).gene_model.tolist()
    bottomgenes=allgenes.filter(matcheslist, list_to_match=downregulated_genes).saveas()
    bottomoperons = alloperons.intersect(bottomgenes, u=True, s=True)
    bottom_operons_upstream=bottomoperons.flank(g=chromosmome_sizes,l=bp_length,r=0,s=True)
    bottom_operons_upstream3=bottom_operons_upstream.filter(lambda x: (x.end-x.start)>100).saveas()
    bottomtranscripts=alltranscripts.filter(matcheslist_mRNA, list_to_match=downregulated_genes).saveas()
    bottomtranscripts_no_operons=bottomtranscripts.subtract(alloperons, s=True, A=True)
    upstream_bottomtranscripts=bottomtranscripts_no_operons.flank(g=chromosmome_sizes,l=bp_length,r=0,s=True)
    upstream_bottomtranscripts3=upstream_bottomtranscripts.filter(lambda x: (x.end-x.start)>100).saveas()
    upstream_bottomtranscripts4=upstream_bottomtranscripts3.sort().merge().saveas()
    bottomenhancers_merged=upstream_bottomtranscripts4.cat(bottom_operons_upstream3, postmerge=True)

    #find final enhancer regions that don't overlap
    final_topenhancers=top_enhancers.subtract(bottomenhancers_merged).saveas()
    final_bottomenhancers=bottomenhancers_merged.subtract(top_enhancers).saveas()

    final_topenhancers.sequence(fi=genome).save_seqs(output_dir + 'target_sequences.fa')
    final_bottomenhancers.sequence(fi=genome).save_seqs(output_dir + 'control_sequences.fa')


    final_top_enhancer_peaks = final_topenhancers.intersect(ATACpeaks, wb=True).saveas()
    final_bottom_enhancer_peaks = final_bottomenhancers.intersect(ATACpeaks, wb=True).saveas()

    out=final_top_enhancer_peaks.sequence(fi=genome).save_seqs(output_dir + 'target_sequences_ATAC_peaks.fa')
    negout=final_bottom_enhancer_peaks.sequence(fi=genome).save_seqs(output_dir + 'control_sequences_ATAC_peaks.fa')

tfbs_enrichment_analysis(sys.argv[1], sys.argv[2], sys.argv[3],
                         sys.argv[4], sys.argv[5], sys.argv[6],
                         sys.argv[7], sys.argv[8], sys.argv[9],
                         sys.argv[10], sys.argv[11])
