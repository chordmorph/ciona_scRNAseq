#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

raw_region_results = pd.read_csv('{}/opossum_output/results.txt'.format(sys.argv[1]), sep = '\t')
ATAC_results = pd.read_csv('{}/opossum_output_ATAC/results.txt'.format(sys.argv[1]), sep = '\t')

joined_data = raw_region_results.merge(ATAC_results, on = 'ID', suffixes=['_raw_regions', '_ATAC_regions'])
joined_data.replace(to_replace = np.inf, value = np.nan, inplace = True)
joined_data = joined_data.assign(average_zscore = joined_data.filter(regex='^Z-score').mean(axis = 1),
                                average_ks_score = joined_data.filter(regex='^KS score').mean(axis = 1),
                                average_fisher_score = joined_data.filter(regex='^Fisher score').mean(axis = 1))
column_order2 = ['ID', 'TF_raw_regions', 'average_zscore', 'Z-score_raw_regions', 'Z-score_ATAC_regions']
joined_data2 = joined_data[column_order2]
joined_data2 = joined_data2.sort_values('average_zscore', ascending=False)
joined_data2 = joined_data2.rename(columns={'TF_raw_regions': 'TF_name'})

joined_data2.to_csv('{}/{}_opossum_average_results.csv'.format(sys.argv[1], sys.argv[2]), index = False)

motif_clusters = pd.read_excel('jaspar_2020_motif_clusters.xlsx')
all_motifs = []
all_clusters_central_motifs = []
all_cluster_ids = []
all_cluster_names = []
num_motifs_in_cluster = []

for i in range(len(motif_clusters.motifs_in_cluster)):
    motif_list = motif_clusters.motifs_in_cluster[i].split(',')
    for motif in motif_list:
        all_motifs.append(motif)
        all_clusters_central_motifs.append(motif_clusters.cluster_motif_id[i])
        all_cluster_ids.append(motif_clusters.cluster_id[i])
        all_cluster_names.append(motif_clusters.cluster_motif_name[i])
        num_motifs_in_cluster.append(motif_clusters.Number_of_Motifs[i])

mapping_df = pd.DataFrame([all_motifs, all_cluster_ids, all_cluster_names, all_clusters_central_motifs, num_motifs_in_cluster]).T
mapping_df.columns = ['TF_name', 'cluster_id', 'cluster_name', 'central_motif_for_cluster', 'cluster_size']

test_data = joined_data2

merged = mapping_df.merge(test_data)
merged = merged.replace(to_replace=np.inf, value=np.nan)
average_zscores = pd.DataFrame(merged.groupby('cluster_id').average_zscore.max())
cluster_names = pd.DataFrame(merged.groupby('cluster_id').cluster_name.max())
cluster_size = pd.DataFrame(merged.groupby('cluster_id').cluster_size.max())
highest_pwm = pd.DataFrame(merged.groupby(['cluster_id']).average_zscore.idxmax())
best_pwm_in_cluster = []

for cluster in range(len(highest_pwm.index)):
    try:
        best_pwm_in_cluster.append(merged.TF_name[int(highest_pwm.average_zscore[cluster])])
    except:
        best_pwm_in_cluster.append('no hits in cluster')
        continue

final_average_results = average_zscores.join([cluster_names, pd.DataFrame(best_pwm_in_cluster, index=highest_pwm.index)])
final_average_results = final_average_results.sort_values(by='average_zscore', ascending=False)
final_average_results = final_average_results.rename(columns= {'average_zscore': 'Max Z-score in Cluster', 'cluster_name': 'Cluster Name',
                                      0: 'Top PWM hit in Cluster'})
final_average_results.to_csv('{}/{}_PWM_cluster_results.csv'.format(sys.argv[1], sys.argv[2]), index = False)
