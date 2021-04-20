#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import scipy.spatial.distance as dist

all_sites_relatedness_data = pd.read_csv('', delimiter = '\t')
all_sites_relatedness_data.RELATEDNESS_AJK.fillna(value=-3.0, inplace=True) #for no overlapping SNPs
all_sites_condensed = all_sites_relatedness_data.RELATEDNESS_AJK[all_sites_relatedness_data.INDV1 != all_sites_relatedness_data.INDV2]
all_sites_relatedness_matrix = dist.squareform(all_sites_condensed)
all_sites_relatedness_linkage = linkage(all_sites_relatedness_matrix, method='ward')


subset_sites_relatedness_data = pd.read_csv('', delimiter = '\t')
subset_sites_relatedness_data.RELATEDNESS_AJK.fillna(value=-3.0, inplace=True) #for no overlapping SNPs
subset_sites_condensed = subset_sites_relatedness_data.RELATEDNESS_AJK[subset_sites_relatedness_data.INDV1 != subset_sites_relatedness_data.INDV2]
subset_sites_relatedness_matrix = dist.squareform(subset_sites_condensed)
subset_sites_relatedness_linkage = linkage(subset_sites_relatedness_matrix, method='ward')

fig1 = plt.figure(figsize=(10,8))
ax1 = fig1.add_subplot(111)
all_sites_mito_relatedness_dend = dendrogram(all_sites_relatedness_linkage, ax = ax1)
fig1.savefig('outputs/all_sites_dendrogram.eps')

fig2 = plt.figure(figsize=(10,8))
ax2 = fig2.add_subplot(111)
subset_sites_mito_relatedness_dend = dendrogram(subset_sites_relatedness_linkage, ax = ax2)
fig2.savefig('outputs/subset_sites_dendrogram.eps')


max_d1 = 350
max_d2 = 250

all_sites_clusters = fcluster(all_sites_relatedness_linkage, max_d1, criterion = 'distance')
subset_sites_clusters = fcluster(subset_sites_relatedness_linkage, max_d2, criterion = 'distance')

barcodes = np.unique(all_sites_relatedness_data.INDV1)
relatedness_output = pd.DataFrame([barcodes, all_sites_clusters, subset_sites_clusters]).T
relatedness_output.columns = ['barcode', 'all_sites_cluster', 'subset_sites_cluster']

relatedness_output.to_csv('', index=False)
